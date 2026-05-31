% =========================================================================
% Monte Carlo com Mapeamento de Envelope
%
% Matheus Martins Marinato — TCC UFSM 2025/2026
%
% SCRIPT quase AUTO-CONTIDO: não precisa rodar run_scurve nem init_mpc antes.   
% RODAR init_ekf ANTES!
%
% DUAS CAMPANHAS:
%
%   CAMPANHA 1 — ENVELOPE:
%     Amostragem ESTRATIFICADA por altitude (200m a 1500m).
%     N_por_faixa runs por faixa, variando as demais CIs.
%     Objetivo: achar onde funciona e onde quebra.
%     Resultado: gráfico de taxa de sucesso vs altitude.
%
%   CAMPANHA 2 — MC DENSO (opcional, desativada por padrão):
%     Estatística dentro do envelope seguro encontrado na Campanha 1.
%     N_denso runs com dispersão realista.
%     Resultado: histogramas, percentis, tabela de requisitos.
%
% USO:
%   Para ativar Campanha 1, mudar RODAR_CAMPANHA_1 = true.
%   Para ativar Campanha 2, mudar RODAR_CAMPANHA_2 = true.
%
% SAÍDA:
%   mc_envelope.mat  — resultados da Campanha 1
%   mc_denso.mat     — resultados da Campanha 2 (se ativada)
%   Ambos contêm TUDO: CIs, estado na ignição, métricas de touchdown.
%
% ESTRUTURA DO ARQUIVO:
%   - Bloco de configuração e pesos (este topo)
%   - Laço da Campanha 1 e da Campanha 2
%   - Funções de infraestrutura (prealoc, geração de CI, pipeline de 1 run)
%   - Funções de relatório e de plots
%   - Funções da S-curve auto-contida
%
% =========================================================================


fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║       MONTE CARLO v2  —  Envelope + Estatística Densa         ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

% =========================================================================
%  CONFIGURAÇÃO GERAL
% =========================================================================

% Modelo Simulink integrado (GNC + planta).
SLX_MODEL = 'Integrated_GNC_and_Rocket_Model_Simulink';

SEED      = 42;
% ── Campanha 1: Envelope ─────────────────────────────────────────────
RODAR_CAMPANHA_1 = true;

% Faixas de altitude amostradas (estratificação) e nº de runs por faixa.
h_faixas    = [200, 300, 400, 500, 600, 800, 1000, 1200, 1500];
N_por_faixa = 15;       % runs por faixa (total = 9×Nporfaixa = ~~ runs)

% ── Campanha 2: MC Denso (dentro do envelope seguro) ─────────────────
RODAR_CAMPANHA_2 = true;

% Limites preenchidos depois de inspecionar o envelope da Campanha 1.
h_denso_min  = 400;          % Preencher após ver resultado da Camp. 1
h_denso_max  = 1200;         % Preencher após ver resultado da Camp. 1
N_denso      = 1000;         % Numero de simulacoes

% ── Modelo de Vento ──────────────────────────────────────────────────
sigma_vento = 5.0;           % [m/s] parâmetro Rayleigh

% =========================================================================
%  PARÂMETROS DO FOGUETE E PESOS
% =========================================================================

p_base = falcon9_params();

%Depois de tuning com bayesopt
% Multiplicadores finais da otimização 2-estágios (mesmos do init_mpc).
k_pos    = 41.8856;
k_vel    = 0.4940;
k_h      = 11.2475;
k_vz     = 1.0737;
k_ang    = 0.4550;
k_rate   = 0.2210;
k_tvc    = 0.2170;

% Matriz de estado Q (Bryson): peso = knob / (desvio_max)^2 por estado.
Q_mpc = diag([...
  k_pos  / (5.0)^2,           ... %  1. x_I
  k_pos  / (5.0)^2,           ... %  2. y_I
  k_h    / (2.0)^2,           ... %  3. z_I
  k_vel  / (3.0)^2,           ... %  4. u_B
  k_vel  / (3.0)^2,           ... %  5. v_B
  k_vz   / (1.0)^2,           ... %  6. w_B
  k_ang  / (deg2rad(5))^2,    ... %  7. phi
  k_ang  / (deg2rad(5))^2,    ... %  8. theta
  k_rate / (deg2rad(5))^2,    ... %  9. p_B
  k_rate / (deg2rad(5))^2     ... % 10. q_B
]);

% Matriz de controle R: empuxo com peso base; TVC ponderado por k_tvc.
R_mpc = diag([...
  1/(p.T_max - p.T_min)^2,    ... % T
  k_tvc / (p.delta_max)^2,     ... % delta_y
  k_tvc / (p.delta_max)^2      ... % delta_z
]);

% Pesos da taxa de variação do comando (rate weighting / anti-chattering).
k_rate_T     = 300;     % multiplicador rate do empuxo
k_rate_delta = 50;    % multiplicador rate do TVC 

S_rate = diag([ ...
    k_rate_T     / (p.T_max - p.T_min)^2, ...    % d(T)/dt
    k_rate_delta / (p.delta_max)^2,        ...    % d(δy)/dt
    k_rate_delta / (p.delta_max)^2          ...    % d(δz)/dt
]);

% Opções do solver QP do MPC (mesmas do init_mpc).
opts_qp = optimoptions('quadprog', ...
    'Display','off', 'Algorithm','interior-point-convex', ...
    'MaxIterations', 50, 'OptimalityTolerance', 1e-6);

% ── Constantes da S-curve ────────────────────────────────────────────
% Parâmetros de propulsão e passo usados pela geração de trajetória.
empuxo_max   = p_base.T_max;
throttle_min = p_base.throttle_min;
T_min_motor  = empuxo_max * throttle_min;
margem       = 0.90;   % usa no máx. 90% do empuxo
dt_sc        = p_base.dt;
g            = p_base.g;
dt = p.dt;

% Condições finais de pouso (fixas)
% Alvo de touchdown: pousar parado e centrado (tudo zero nos três canais).
hf=0; vzf=0; azf=0; xf=0; vxf=0; axf=0; yf=0; vyf=0; ayf=0;

% Parâmetros do scan
% SCAN_DT: passo do scan de ignição; TVIOL: tolerância de violação de T_min;
% ALPHA: posição relativa do instante de ignição dentro da janela viável.
SCAN_DT  = 0.2;
TVIOL    = 0.03;
ALPHA    = 0.5;

% ── Perturbações das CIs (tudo EXCETO altitude, que é estratificada) ─
% Desvios-padrão (1σ) aplicados sobre os valores nominais em cada run.
pert.vz0_sigma    = 5.0;             % [m/s]
pert.xy_sigma     = 8.0;             % [m]
pert.phi_sigma    = deg2rad(2.0);    % [rad]
pert.theta_sigma  = deg2rad(2.0);    % [rad]
pert.q0_sigma     = deg2rad(0.5);    % [rad/s]
pert.m0_sigma     = 500;             % [kg]

% Valores nominais (exceto h0 que vem da faixa)
nom.vz0    = -50.0;
nom.x0     = 10.0;
nom.y0     = -8.0;
nom.phi0   = deg2rad(-3.0);
nom.theta0 = deg2rad(3.0);
nom.q0     = deg2rad(0.5);
nom.m0     = p_base.m0;

% ── Carregar Simulink ────────────────────────────────────────────────
% Carrega o modelo uma vez; cada run reaproveita a mesma sessão carregada.
fprintf('Carregando modelo Simulink...');
load_system(SLX_MODEL);
fprintf(' OK\n\n');

% #########################################################################
%  CAMPANHA 1 — MAPEAMENTO DE ENVELOPE
% #########################################################################

if RODAR_CAMPANHA_1

N_faixas = length(h_faixas);
N_total  = N_faixas * N_por_faixa;

fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  CAMPANHA 1 — Mapeamento de Envelope\n');
fprintf('  Faixas: %s m\n', mat2str(h_faixas));
fprintf('  Runs/faixa: %d  |  Total: %d runs\n', N_por_faixa, N_total);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

rng(SEED);  

% ── Pré-alocar struct de resultados ──────────────────────────────────
R = prealoc_resultados(N_total);

fprintf('  Run  | h0_nom |  Status  |  err_lat  |  vz_td  | vlat_td | theta_max | Tempo  |  ETA\n');
fprintf('  -----+--------+----------+-----------+---------+---------+-----------+--------+-------\n');

t_mc = tic;
run_idx = 0;

% Laço externo nas faixas de altitude, interno nas repetições por faixa.
for fi = 1:N_faixas
    h0_nom = h_faixas(fi);

    for ri = 1:N_por_faixa
        run_idx = run_idx + 1;
        t_run = tic;

        % ── Gerar CIs ───────────────────────────────────────────────
        % Altitude: uniforme dentro da faixa ±15% pra dar variação
        % (mantém a estratificação, mas evita repetir exatamente o mesmo h0).
        h0_spread = h0_nom * 0.15;
        h0_k = h0_nom + h0_spread * (2*rand()-1);
        h0_k = max(50, h0_k);

        % Sorteia as demais CIs e o vento em torno dos valores nominais.
        [ic, vw] = gerar_ci(h0_k, nom, pert, sigma_vento);

        % Salvar CIs
        R = salvar_ic(R, run_idx, ic, vw, h0_nom);

        % ── Rodar uma simulação ─────────────────────────────────────
        % Pipeline completo: gera S-curve, monta x0, simula e extrai métricas.
        [status, metricas, msg, ign] = rodar_uma_sim( ...
            ic, vw, p_base, Q_lqr, R_lqr, opts_qp, ...
            empuxo_max, T_min_motor, margem, dt_sc, g, ...
            SCAN_DT, TVIOL, ALPHA, ...
            hf,vzf,azf, xf,vxf,axf, yf,vyf,ayf, ...
            SLX_MODEL);

        % Salvar resultados
        R.status(run_idx)    = status;
        R.msg{run_idx}       = msg;
        R.faixa_idx(run_idx) = fi;
        R.h0_nominal(run_idx)= h0_nom;

        % Registra o estado na ignição (se a S-curve foi gerada).
        if ~isnan(ign.h_ign)
            R.ign.h(run_idx)     = ign.h_ign;
            R.ign.vz(run_idx)    = ign.vz_ign;
            R.ign.x(run_idx)     = ign.x_ign;
            R.ign.y(run_idx)     = ign.y_ign;
            R.ign.T_burn(run_idx)= ign.T_burn;
            R.ign.t_ign(run_idx) = ign.t_ign;
        end

        if status >= 0  % sim rodou (sucesso ou crash)
            R = salvar_metricas(R, run_idx, metricas);
        end

        % ── Progresso ──────────────────────────────────────────────
        % Imprime a linha da run com status, métricas e ETA.
        dt_run = toc(t_run);
        vlat_k = sqrt(R.vx_td(run_idx)^2 + R.vy_td(run_idx)^2);
        print_progress_v2(run_idx, N_total, h0_nom, status, ...
            R.err_lat(run_idx), R.vz_td(run_idx), vlat_k, ...
            R.theta_max(run_idx), dt_run, t_mc);
    end
end

t_total_c1 = toc(t_mc);

% ── Salvar ───────────────────────────────────────────────────────────
% Empacota resultados + configuração da campanha e grava em .mat (-v7.3
% para suportar structs grandes).
campanha1.R          = R;
campanha1.h_faixas   = h_faixas;
campanha1.N_por_faixa= N_por_faixa;
campanha1.pert       = pert;
campanha1.nom        = nom;
campanha1.sigma_vento= sigma_vento;
campanha1.SEED       = SEED;
campanha1.tempo_min  = t_total_c1 / 60;

save('mc_envelope.mat', 'campanha1', '-v7.3');
fprintf('\n  Resultados salvos em: mc_envelope.mat\n');

% ── Relatório do Envelope ────────────────────────────────────────────
imprimir_relatorio_envelope(R, h_faixas, N_por_faixa, t_total_c1);

% ── Plots do Envelope ────────────────────────────────────────────────
plotar_envelope(R, h_faixas, N_por_faixa, sigma_vento);
salvar_crashes(R,  'crashes_envelope.csv');
end  % CAMPANHA 1


% #########################################################################
%  CAMPANHA 2 — MC DENSO DENTRO DO ENVELOPE
% #########################################################################

if RODAR_CAMPANHA_2

fprintf('\n\n');
fprintf('═══════════════════════════════════════════════════════════════\n');
fprintf('  CAMPANHA 2 — MC Denso [%d, %d] m\n', h_denso_min, h_denso_max);
fprintf('  Runs: %d\n', N_denso);
fprintf('═══════════════════════════════════════════════════════════════\n\n');

rng(SEED + 1000);   % semente diferente da Camp. 1

R2 = prealoc_resultados(N_denso);
fprintf('  Run  | h0_nom |  Status  |  err_lat  |  vz_td  | vlat_td | theta_max | Tempo  |  ETA\n');
fprintf('  -----+--------+----------+-----------+---------+---------+-----------+--------+-------\n');

t_mc2 = tic;

% Sem estratificação aqui: altitude uniforme em todo o envelope seguro.
for k = 1:N_denso
    t_run = tic;

    % Altitude uniforme dentro do envelope seguro
    h0_k = h_denso_min + (h_denso_max - h_denso_min) * rand();

    [ic, vw] = gerar_ci(h0_k, nom, pert, sigma_vento);
    R2 = salvar_ic(R2, k, ic, vw, h0_k);

    [status, metricas, msg, ign] = rodar_uma_sim( ...
        ic, vw, p_base, Q_lqr, R_lqr, opts_qp, ...
        empuxo_max, T_min_motor, margem, dt_sc, g, ...
        SCAN_DT, TVIOL, ALPHA, ...
        hf,vzf,azf, xf,vxf,axf, yf,vyf,ayf, ...
        SLX_MODEL);

    R2.status(k) = status;
    R2.msg{k}    = msg;

    if ~isnan(ign.h_ign)
        R2.ign.h(k)     = ign.h_ign;
        R2.ign.vz(k)    = ign.vz_ign;
        R2.ign.x(k)     = ign.x_ign;
        R2.ign.y(k)     = ign.y_ign;
        R2.ign.T_burn(k)= ign.T_burn;
        R2.ign.t_ign(k) = ign.t_ign;
    end

    if status >= 0
        R2 = salvar_metricas(R2, k, metricas);
    end

    dt_run = toc(t_run);
    vlat_k = sqrt(R2.vx_td(k)^2 + R2.vy_td(k)^2);
    print_progress_v2(k, N_denso, h0_k, status, ...
        R2.err_lat(k), R2.vz_td(k), vlat_k, R2.theta_max(k), dt_run, t_mc2);
end

t_total_c2 = toc(t_mc2);

% ── Salvar ───────────────────────────────────────────────────────────
campanha2.R           = R2;
campanha2.h_range     = [h_denso_min, h_denso_max];
campanha2.N           = N_denso;
campanha2.pert        = pert;
campanha2.nom         = nom;
campanha2.sigma_vento = sigma_vento;
campanha2.SEED        = SEED + 1000;
campanha2.tempo_min   = t_total_c2 / 60;

save('mc_denso.mat', 'campanha2', '-v7.3');
fprintf('\n  Resultados salvos em: mc_denso.mat\n');

% ── Relatório e Plots ────────────────────────────────────────────────
imprimir_relatorio_denso(R2, h_denso_min, h_denso_max, t_total_c2);
plotar_mc_denso(R2, h_denso_min, h_denso_max, sigma_vento);
salvar_crashes(R2, 'crashes_dense.csv'); 
end  % CAMPANHA 2

fprintf('\n  Monte Carlo v2 concluído.\n\n');


% #########################################################################
% #########################################################################
%  FUNÇÕES: INFRAESTRUTURA DO MC
% #########################################################################
% #########################################################################

% =========================================================================
%  prealoc_resultados — cria struct vazia pra N runs
% =========================================================================
function R = prealoc_resultados(N)
% Aloca todos os campos de resultado com NaN/zeros para N runs. Pré-alocar
% evita realocação dentro do laço e deixa claro o esquema de dados coletados.
    R.status     = zeros(N,1);     % 1=ok, 0=crash, -1=scurve fail, -2=sim fail
    R.msg        = cell(N,1);
    R.faixa_idx  = zeros(N,1);
    R.h0_nominal = NaN(N,1);

    % Condições iniciais
    R.ic.h0     = NaN(N,1);
    R.ic.vz0    = NaN(N,1);
    R.ic.x0     = NaN(N,1);
    R.ic.y0     = NaN(N,1);
    R.ic.phi0   = NaN(N,1);
    R.ic.theta0 = NaN(N,1);
    R.ic.q0     = NaN(N,1);
    R.ic.m0     = NaN(N,1);
    R.ic.vw_x   = NaN(N,1);
    R.ic.vw_y   = NaN(N,1);

    % Estado na ignição
    R.ign.h      = NaN(N,1);
    R.ign.vz     = NaN(N,1);
    R.ign.x      = NaN(N,1);
    R.ign.y      = NaN(N,1);
    R.ign.T_burn = NaN(N,1);
    R.ign.t_ign  = NaN(N,1);

    % Métricas de touchdown
    R.x_td      = NaN(N,1);
    R.y_td      = NaN(N,1);
    R.z_td      = NaN(N,1);
    R.vx_td     = NaN(N,1);
    R.vy_td     = NaN(N,1);
    R.vz_td     = NaN(N,1);
    R.phi_td    = NaN(N,1);
    R.theta_td  = NaN(N,1);
    R.psi_td    = NaN(N,1);
    R.err_lat   = NaN(N,1);
    R.t_sim     = NaN(N,1);
    R.theta_max = NaN(N,1);
    R.phi_max   = NaN(N,1);
end

% =========================================================================
%  gerar_ci — gera condições iniciais aleatórias
% =========================================================================
function [ic, vw] = gerar_ci(h0, nom, pert, sigma_vento)
% Sorteia uma condição inicial: cada estado = nominal + ruído gaussiano
% truncado (tg). h0 vem de fora (estratificado). Devolve também o vento.
    ic.h0     = h0;
    ic.vz0    = nom.vz0    + tg(pert.vz0_sigma);
    ic.x0     = nom.x0     + tg(pert.xy_sigma);
    ic.y0     = nom.y0     + tg(pert.xy_sigma);
    ic.phi0   = nom.phi0   + tg(pert.phi_sigma);
    ic.theta0 = nom.theta0 + tg(pert.theta_sigma);
    ic.q0     = nom.q0     + tg(pert.q0_sigma);
    ic.m0     = nom.m0     + tg(pert.m0_sigma);

    % Sanidade
    % Garante descida (vz0 negativo) e massa dentro de limites físicos.
    ic.vz0 = min(-10, ic.vz0);
    ic.m0  = max(22500 + 2000, min(30000 + 3000, ic.m0));

    % Vento Rayleigh + direção uniforme
    % Magnitude com distribuição de Rayleigh; direção uniforme em [0,2π).
    vw_mag = sigma_vento * sqrt(-2*log(rand()));
    vw_dir = 2*pi*rand();
    vw.x   = vw_mag * cos(vw_dir);
    vw.y   = vw_mag * sin(vw_dir);
end

% =========================================================================
%  salvar_ic — registra CIs na struct
% =========================================================================
function R = salvar_ic(R, k, ic, vw, h0_nom)
% Copia a CI sorteada (e o vento) para a linha k da struct de resultados.
    R.ic.h0(k)     = ic.h0;
    R.ic.vz0(k)    = ic.vz0;
    R.ic.x0(k)     = ic.x0;
    R.ic.y0(k)     = ic.y0;
    R.ic.phi0(k)   = ic.phi0;
    R.ic.theta0(k) = ic.theta0;
    R.ic.q0(k)     = ic.q0;
    R.ic.m0(k)     = ic.m0;
    R.ic.vw_x(k)   = vw.x;
    R.ic.vw_y(k)   = vw.y;
    R.h0_nominal(k)= h0_nom;
end

% =========================================================================
%  salvar_metricas — registra métricas de touchdown na struct
% =========================================================================
function R = salvar_metricas(R, k, m)
% Copia as métricas de touchdown (struct m) para a linha k da struct R.
    R.x_td(k)      = m.x_td;
    R.y_td(k)      = m.y_td;
    R.z_td(k)      = m.z_td;
    R.vx_td(k)     = m.vx_td;
    R.vy_td(k)     = m.vy_td;
    R.vz_td(k)     = m.vz_td;
    R.phi_td(k)    = m.phi_td;
    R.theta_td(k)  = m.theta_td;
    R.psi_td(k)    = m.psi_td;
    R.err_lat(k)   = m.err_lat;
    R.t_sim(k)     = m.t_sim;
    R.theta_max(k) = m.theta_max;
    R.phi_max(k)   = m.phi_max;
end

% =========================================================================
%  rodar_uma_sim — pipeline completo de uma run
% =========================================================================
function [status, metricas, msg, ign] = rodar_uma_sim( ...
    ic, vw, p_base, Q_lqr, R_lqr, opts_qp, ...
    empuxo_max, T_min_motor, margem, dt_sc, g, ...
    SCAN_DT, TVIOL, ALPHA, ...
    hf,vzf,azf, xf,vxf,axf, yf,vyf,ayf, ...
    SLX_MODEL)
% Executa uma run completa e devolve o desfecho codificado em 'status':
%   1 = pouso bem-sucedido, 0 = crash (chegou ao solo violando requisito),
%  -1 = falha ao gerar a S-curve, -2 = falha na simulação/extração.

    metricas = struct();
    msg = '';
    ign = struct('h_ign',NaN,'vz_ign',NaN,'x_ign',NaN,'y_ign',NaN,...
                 'T_burn',NaN,'t_ign',NaN);

    % ── 1. Gerar S-curve ────────────────────────────────────────────
    % Falha aqui (envelope fechado, sem impacto) → status -1, aborta a run.
    try
        [traj, ign_out] = gerar_scurve_v2( ...
            ic.h0, ic.vz0, ic.x0, ic.y0, ic.phi0, ic.theta0, ...
            hf,vzf,azf, xf,vxf,axf, yf,vyf,ayf, ...
            empuxo_max, T_min_motor, ic.m0, margem, dt_sc, g, ...
            SCAN_DT, TVIOL, ALPHA);
        ign = ign_out;
    catch ME
        status = -1;
        msg = ME.message;
        return;
    end

    % ── 2. Montar x0_sim ────────────────────────────────────────────
    % Estado inicial de 13 componentes a partir da CI sorteada.
    x0_sim = [ic.x0; ic.y0; ic.h0; 0; 0; ic.vz0; ...
              ic.phi0; ic.theta0; 0; 0; ic.q0; 0; ic.m0];

    % ── 3. Colocar no workspace ─────────────────────────────────────
    % O modelo Simulink lê estas variáveis do base workspace via SrcWorkspace.
    assignin('base', 'p',        p_base);
    assignin('base', 'x0_sim',   x0_sim);
    assignin('base', 'Q_lqr',    Q_lqr);
    assignin('base', 'R_lqr',    R_lqr);
    assignin('base', 'opts_qp',  opts_qp);
    assignin('base', 'traj',     traj);
    assignin('base', 'v_vento_mc', [vw.x; vw.y; 0]);

    % ── 4. Simular ──────────────────────────────────────────────────
    % Falha de simulação → status -2.
    try
        sim_out = sim(SLX_MODEL, 'SrcWorkspace', 'base', 'StopTime', '60');
    catch ME
        status = -2;
        msg = ['Sim falhou: ' ME.message];
        return;
    end

    % ── 5. Extrair métricas ─────────────────────────────────────────
    % x_sim pode voltar como timeseries, struct-with-time ou matriz; o bloco
    % abaixo normaliza para [amostras × 13] antes de ler o touchdown.
    try
        ts = sim_out.get('x_sim');
        if isa(ts, 'timeseries')
            t_vec = ts.Time;
            dados = squeeze(ts.Data);
            if size(dados,1) ~= length(t_vec), dados = dados'; end
        elseif isstruct(ts)
            t_vec = ts.time;
            dados = ts.signals.values;
        else
            dados = double(ts);
            t_vec = sim_out.get('tout');
            if isa(t_vec,'timeseries'), t_vec = t_vec.Data; end
        end

        % Estado final = touchdown; base de todas as métricas.
        xf_k  = dados(end,:);
        t_end = t_vec(end);

        metricas.x_td     = xf_k(1);
        metricas.y_td     = xf_k(2);
        metricas.z_td     = xf_k(3);
        metricas.vx_td    = xf_k(4);
        metricas.vy_td    = xf_k(5);
        metricas.vz_td    = xf_k(6);
        metricas.phi_td   = rad2deg(xf_k(7));
        metricas.theta_td = rad2deg(xf_k(8));
        metricas.psi_td   = rad2deg(xf_k(9));
        metricas.err_lat  = sqrt(xf_k(1)^2 + xf_k(2)^2);
        metricas.t_sim    = t_end;
        metricas.theta_max= max(abs(rad2deg(dados(:,8))));
        metricas.phi_max  = max(abs(rad2deg(dados(:,7))));

        % ── Success criterion PFC eq. (8.1) ────────────────────────────
    %   |w_td|   <= 7.50 m/s   (structural sink-rate limit)
    %   ||r_td|| <= 10 m      (geometric containment within landing zone)
    % Plus z < 5 m as a sanity check that the run actually reached the pad
    
    
    % status 1 = pouso válido; 0 = chegou ao solo mas violou algum critério.
    if xf_k(3) < 5.0 && abs(xf_k(6)) <= 7.5 && metricas.err_lat <= 10.0
        status = 1;
    else
        status = 0;
        msg = sprintf('fail: z=%.1f w=%.1f rlat=%.1f', ...
                      xf_k(3), xf_k(6), metricas.err_lat);
    end

    catch ME
        status = -2;
        msg = ['Extração falhou: ' ME.message];
    end
end


% #########################################################################
%  FUNÇÕES: RELATÓRIOS
% #########################################################################

% =========================================================================
%  imprimir_relatorio_envelope
% =========================================================================
function imprimir_relatorio_envelope(R, h_faixas, N_por_faixa, t_total)
% Imprime, faixa a faixa, a taxa de sucesso e as métricas medianas, e ao
% final sugere o envelope seguro (faixas com taxa >= 90%) para a Campanha 2.
    N_faixas = length(h_faixas);
    N_total  = length(R.status);

    fprintf('\n\n');
    fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
    fprintf('║              RELATÓRIO — CAMPANHA 1 (ENVELOPE)                ║\n');
    fprintf('╠═══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Runs total: %d   |  Tempo: %.1f min                       \n', N_total, t_total/60);
    fprintf('║  Sucessos:   %d (%.1f%%)  |  S-fail: %d  |  Crash: %d      \n', ...
            sum(R.status==1), 100*sum(R.status==1)/N_total, ...
            sum(R.status==-1), sum(R.status==0));
    fprintf('╠═══════════════════════════════════════════════════════════════╣\n');
    fprintf('║                                                               ║\n');
    fprintf('║  ENVELOPE DE OPERAÇÃO POR FAIXA DE ALTITUDE                   ║\n');
    fprintf('║                                                               ║\n');
    fprintf('║  h0 [m]  | N_ok/N | Taxa  | err_lat_med | |vz|_med | Nota     ║\n');
    fprintf('║  --------+--------+-------+-------------+----------+-------   ║\n');

    % Para cada faixa: conta sucessos e calcula medianas só sobre os sucessos.
    for fi = 1:N_faixas
        mask = (R.faixa_idx == fi);
        n_total_f = sum(mask);
        n_ok_f    = sum(mask & R.status == 1);
        taxa      = 100 * n_ok_f / max(1, n_total_f);

        if n_ok_f > 0
            ok_mask_f = mask & R.status == 1;
            el_med    = median(R.err_lat(ok_mask_f));
            vz_med    = median(abs(R.vz_td(ok_mask_f)));
        else
            el_med = NaN;
            vz_med = NaN;
        end

        % Classificação qualitativa da faixa pela taxa de sucesso.
        if taxa >= 90
            nota = ' BOM  ';
        elseif taxa >= 50
            nota = 'MARGIN';
        elseif taxa > 0
            nota = 'FRACO ';
        else
            nota = 'IMPOSSIVEL';
        end

        fprintf('║  %6d  | %2d/%2d  | %5.1f%% | %9.2f m | %6.2f   | %s\n', ...
                h_faixas(fi), n_ok_f, n_total_f, taxa, el_med, vz_med, nota);
    end

    fprintf('║                                                              ║\n');

    % ── Sugestão de envelope ─────────────────────────────────────────
    % Recalcula as taxas e identifica a faixa contígua com taxa >= 90%.
    taxas = zeros(N_faixas,1);
    for fi = 1:N_faixas
        mask = (R.faixa_idx == fi);
        n_t = sum(mask);
        n_o = sum(mask & R.status == 1);
        taxas(fi) = n_o / max(1, n_t);
    end

    idx_bom = find(taxas >= 0.9);
    if ~isempty(idx_bom)
        h_min_env = h_faixas(idx_bom(1));
        h_max_env = h_faixas(idx_bom(end));
        fprintf('║  ENVELOPE SEGURO (taxa >= 90%%): [%d, %d] m               \n', h_min_env, h_max_env);
        fprintf('║  → Use estes valores em h_denso_min/max da Campanha 2.    ║\n');
    else
        fprintf('║  [AVISO] Nenhuma faixa atingiu 90%% de sucesso.           ║\n');
    end

    fprintf('╚═══════════════════════════════════════════════════════════════╝\n');
end

% =========================================================================
%  imprimir_relatorio_denso
% =========================================================================
function imprimir_relatorio_denso(R, h_min, h_max, t_total)
% Estatística agregada da Campanha 2: medianas/p95/máximos das métricas de
% touchdown e o confronto p95-vs-requisito (tabela de verificação).
    N = length(R.status);
    ok = (R.status == 1);
    n_ok = sum(ok);

    fprintf('\n\n');
    fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
    fprintf('║              RELATÓRIO — CAMPANHA 2 (MC DENSO)                ║\n');
    fprintf('╠═══════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Runs: %d  |  Envelope: [%d, %d] m  |  Tempo: %.1f min     \n', N, h_min, h_max, t_total/60);
    fprintf('║  Sucessos: %d (%.1f%%)  |  S-fail: %d  |  Crash: %d        \n', ...
            n_ok, 100*n_ok/N, sum(R.status==-1), sum(R.status==0));
    fprintf('╠═══════════════════════════════════════════════════════════════╣\n');

    % Só agrega se houver amostra mínima de sucessos.
    if n_ok > 2
        % Vetores das métricas sobre os sucessos.
        el   = R.err_lat(ok);
        vz   = abs(R.vz_td(ok));
        vlat = sqrt(R.vx_td(ok).^2 + R.vy_td(ok).^2);
        att  = sqrt(R.phi_td(ok).^2 + R.theta_td(ok).^2);
        tm   = R.theta_max(ok);

        fprintf('║  Erro lateral:     med=%.2f  p95=%.2f  max=%.2f m       \n', median(el), prctile(el,95), max(el));
        fprintf('║  |w_B| touchdown:  med=%.2f  p95=%.2f  max=%.2f m/s     \n', median(vz), prctile(vz,95), max(vz));
        fprintf('║  v_lat touchdown:  med=%.2f  p95=%.2f  max=%.2f m/s     \n', median(vlat), prctile(vlat,95), max(vlat));
        fprintf('║  |att| touchdown:  med=%.2f  p95=%.2f  max=%.2f°        \n', median(att), prctile(att,95), max(att));
        fprintf('║  theta_max voo:    med=%.2f  p95=%.2f  max=%.2f°        \n', median(tm), prctile(tm,95), max(tm));
        fprintf('╠═══════════════════════════════════════════════════════════╣\n');

        % Confronta o p95 de cada métrica com seu requisito (chk → OK/FALHA).
        fprintf('║  VERIFICAÇÃO DE REQUISITOS:                              ║\n');
        fprintf('║  %-22s  %7s  %7s  %s\n', 'Métrica', 'p95', 'Req', 'Status');
        chk('Erro lateral [m]',     prctile(el,95),   5.0);
        chk('|vz| td [m/s]',        prctile(vz,95),   4.0);
        chk('Vel. lateral [m/s]',   prctile(vlat,95), 1.5); 
        chk('|att| td [deg]',       prctile(att,95),  6.0);
        chk('theta_max [deg]',      prctile(tm,95),  25.0);
        fprintf('║  %-22s  %6.1f%%  %6.1f%%  %s\n', 'Taxa sucesso', ...
                100*n_ok/N, 90.0, iff(100*n_ok/N>=90,'  OK','  FALHA'));
    end

    fprintf('╚═══════════════════════════════════════════════════════════════╝\n');
end


% #########################################################################
%  FUNCTIONS: PLOTS
% #########################################################################
%
% =========================================================================
%  plotar_envelope — Campaign 1 plots
% =========================================================================
function plotar_envelope(R, h_faixas, N_por_faixa, sigma_vento)
% Painel 2x3 da Campanha 1: taxa de sucesso, erro e velocidade por faixa,
% dispersão dos touchdowns e composição de desfechos por faixa.
    N_faixas = length(h_faixas);
    
    % Compute per-band statistics
    % Estatísticas por faixa (taxa de sucesso, mediana e p95 de erro e vz).
    taxa_suc  = zeros(N_faixas,1);
    el_med    = NaN(N_faixas,1);
    el_p95    = NaN(N_faixas,1);
    vz_med    = NaN(N_faixas,1);
    vz_p95    = NaN(N_faixas,1);
    
    for fi = 1:N_faixas
        mask  = (R.faixa_idx == fi);
        n_t   = sum(mask);
        n_ok  = sum(mask & R.status == 1);
        taxa_suc(fi) = 100 * n_ok / max(1, n_t);
        
        ok_f = mask & R.status == 1;
        if sum(ok_f) > 0
            el_med(fi) = median(R.err_lat(ok_f));
            el_p95(fi) = prctile(R.err_lat(ok_f), 95);
            vz_med(fi) = median(abs(R.vz_td(ok_f)));
            vz_p95(fi) = prctile(abs(R.vz_td(ok_f)), 95);
        end
    end
    
    figure('Name','Campaign 1 — Operational Envelope', ...
           'NumberTitle','off','Position',[50 50 1400 750]);
           
    % ── 1. Success rate vs altitude ──────────────────────────────────
    % Barras de taxa de sucesso; linha de 90% marca o critério de envelope.
    subplot(2,3,1);
    bar(h_faixas, taxa_suc, 0.6, 'FaceColor', [0.3 0.6 0.9]); hold on;
    yline(90, 'r--', '90%', 'LineWidth', 1.5);
    ylim([0 110]);
    xlabel('Nominal h_0 [m]'); ylabel('Success Rate [%]');
    title('Envelope: Success Rate'); grid on;
    
    % Color bars by performance band
    % Recolore cada barra conforme a faixa de desempenho (verde/amarelo/vermelho).
    for fi = 1:N_faixas
        if taxa_suc(fi) >= 90
            % already blue
        elseif taxa_suc(fi) >= 50
            bar(h_faixas(fi), taxa_suc(fi), 0.6, 'FaceColor', [0.9 0.7 0.2]);
        else
            bar(h_faixas(fi), taxa_suc(fi), 0.6, 'FaceColor', [0.9 0.3 0.3]);
        end
    end
    
    % ── 2. Lateral error per band ────────────────────────────────────
    % Erro lateral: mediana e p95 por faixa.
    subplot(2,3,2);
    plot(h_faixas, el_med, 'bo-', 'LineWidth', 1.8, 'MarkerSize', 8); hold on;
    plot(h_faixas, el_p95, 'r^--', 'LineWidth', 1.5, 'MarkerSize', 7);
    xlabel('h_0 [m]'); ylabel('Lateral error [m]');
    title('Lateral Error vs Altitude');
    legend('Median', 'p95', 'Location', 'best'); grid on;
    
    % ── 3. Touchdown velocity per band ───────────────────────────────
    % Velocidade vertical de toque: mediana e p95 por faixa.
    subplot(2,3,3);
    plot(h_faixas, vz_med, 'bo-', 'LineWidth', 1.8, 'MarkerSize', 8); hold on;
    plot(h_faixas, vz_p95, 'r^--', 'LineWidth', 1.5, 'MarkerSize', 7);
    xlabel('h_0 [m]'); ylabel('Touchdown |w_B| [m/s]');
    title('Touchdown Velocity vs Altitude');
    legend('Median', 'p95', 'Location', 'best'); grid on;
    
    % ── 4. Scatter of all touchdowns (color = status) ────────────────
    % Dispersão dos pontos de toque; cor = altitude inicial, X = crash.
    % Círculos verde/vermelho marcam as zonas de 2 m e 5 m em torno do alvo.
    subplot(2,3,4);
    ok   = R.status == 1;
    fail = R.status == 0;
    sfail= R.status == -1;
    
    if any(ok)
        scatter(R.x_td(ok), R.y_td(ok), 30, R.ic.h0(ok), 'filled', ...
                'MarkerFaceAlpha', 0.7); hold on;
    end
    if any(fail)
        scatter(R.x_td(fail), R.y_td(fail), 60, 'r', 'x', 'LineWidth',2); hold on;
    end
    
    th = linspace(0,2*pi,100);
    plot(2*cos(th),2*sin(th),'g--','LineWidth',1.5);
    plot(5*cos(th),5*sin(th),'r--','LineWidth',1.5);
    plot(0,0,'r+','MarkerSize',15,'LineWidth',2);
    
    colorbar; colormap(gca,'parula');
    axis equal; grid on;
    xlabel('X [m]'); ylabel('Y [m]');
    title('Touchdown (color = h_0)');
    
    % ── 5. h0 vs err_lat with color = |vz| ───────────────────────────
    % Erro lateral vs altitude real, colorido pela velocidade de toque.
    subplot(2,3,5);
    if any(ok)
        scatter(R.ic.h0(ok), R.err_lat(ok), 40, abs(R.vz_td(ok)), 'filled'); hold on;
        colorbar; colormap(gca,'hot');
    end
    if any(sfail)
        scatter(R.ic.h0(sfail), zeros(sum(sfail),1), 50, 'r', 'x', 'LineWidth', 2);
    end
    xlabel('Actual h_0 [m]'); ylabel('Lateral error [m]');
    title('h_0 vs Error (color = |vz_{td}|)'); grid on;
    
    % ── 6. Failure type per band (stacked bar) ───────────────────────
    % Barras empilhadas com a composição de desfechos (sucesso/crash/falhas).
    subplot(2,3,6);
    cnt = zeros(N_faixas, 4);   % [success, crash, scurve_fail, sim_fail]
    for fi = 1:N_faixas
        mask = (R.faixa_idx == fi);
        cnt(fi,1) = sum(mask & R.status == 1);
        cnt(fi,2) = sum(mask & R.status == 0);
        cnt(fi,3) = sum(mask & R.status == -1);
        cnt(fi,4) = sum(mask & R.status == -2);
    end
    
    b = bar(h_faixas, cnt, 'stacked');
    b(1).FaceColor = [0.3 0.7 0.4];
    b(2).FaceColor = [0.9 0.3 0.3];
    b(3).FaceColor = [0.9 0.7 0.2];
    b(4).FaceColor = [0.5 0.5 0.5];
    
    xlabel('Nominal h_0 [m]'); ylabel('Count');
    title('Outcome per Band');
    legend('Success','Crash','S-curve fail','Sim fail','Location','best');
    grid on;
    
    sgtitle(sprintf('Campaign 1 — Envelope  |  %d runs  |  wind \\sigma=%.0f m/s', ...
            length(R.status), sigma_vento), 'FontSize', 12, 'FontWeight', 'bold');
end

% =========================================================================
%  plotar_mc_denso — Campaign 2 plots (TWO FIGURES)
% =========================================================================
function plotar_mc_denso(R, h_min, h_max, sigma_vento)
% Duas figuras da Campanha 2: (1) métricas clássicas de touchdown
% (dispersão, histogramas, boxplot) e (2) dinâmica/esforço de controle.
    ok = R.status == 1;
    if sum(ok) < 3, return; end   % sem amostra suficiente, não plota

    % Métricas sobre os sucessos.
    el   = R.err_lat(ok);
    vz   = abs(R.vz_td(ok));
    vlat = sqrt(R.vx_td(ok).^2 + R.vy_td(ok).^2);
    att  = sqrt(R.phi_td(ok).^2 + R.theta_td(ok).^2);

    % ─────────────────────────────────────────────────────────────────────
    % FIGURE 1: Touchdown Metrics (Classical)
    % ─────────────────────────────────────────────────────────────────────
    figure('Name','Campaign 2 — Dense MC (Touchdown)','NumberTitle','off',...
           'Position',[60 60 1400 700]);

    % Dispersão dos pontos de toque com as zonas de 2 m e 5 m.
    subplot(2,3,1);
    scatter(R.x_td(ok), R.y_td(ok), 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5); hold on;
    th = linspace(0,2*pi,100);
    plot(2*cos(th),2*sin(th),'g--','LineWidth',1.5);
    plot(5*cos(th),5*sin(th),'r--','LineWidth',1.5);
    plot(0,0,'r+','MarkerSize',15,'LineWidth',2);
    fail = R.status == 0;
    if any(fail)
        scatter(R.x_td(fail), R.y_td(fail), 60, 'r', 'x', 'LineWidth', 2);
    end
    axis equal; grid on;
    xlabel('X [m]'); ylabel('Y [m]'); title('Touchdown Position');

    % Histograma do erro lateral, com p95 e requisito.
    subplot(2,3,2);
    histogram(el, 20, 'FaceColor', [0.3 0.6 0.9]); hold on;
    xline(prctile(el,95), 'r--', sprintf('p95=%.1fm',prctile(el,95)), 'LineWidth', 1.5);
    xline(5, 'k--', 'Req 5m'); xlabel('Lateral Error [m]'); ylabel('Count');
    title('Lateral Error'); grid on;

    % Histograma da velocidade vertical de toque.
    subplot(2,3,3);
    histogram(vz, 20, 'FaceColor', [0.9 0.5 0.3]); hold on;
    xline(prctile(vz,95), 'r--', sprintf('p95=%.1f',prctile(vz,95)), 'LineWidth', 1.5);
    xline(4, 'k--', 'Req'); xlabel('|w_B| [m/s]'); ylabel('Count');
    title('Vertical Velocity'); grid on;

    % Histograma da atitude de toque.
    subplot(2,3,4);
    histogram(att, 20, 'FaceColor', [0.5 0.8 0.5]); hold on;
    xline(prctile(att,95), 'r--', sprintf('p95=%.1f°',prctile(att,95)), 'LineWidth', 1.5);
    xline(6, 'k--', 'Req'); xlabel('|Attitude| [deg]'); ylabel('Count');
    title('Touchdown Attitude'); grid on;

    % Erro lateral vs altitude, colorido pela velocidade vertical.
    subplot(2,3,5);
    scatter(R.ic.h0(ok), R.err_lat(ok), 30, abs(R.vz_td(ok)), 'filled');
    colorbar; colormap(gca,'hot');
    xlabel('h_0 [m]'); ylabel('Lateral error [m]');
    title('h_0 vs Error (color = |vz|)'); grid on;

    % Boxplot resumo da dispersão das três métricas principais.
    subplot(2,3,6);
    data_box = [el, vz, att];
    boxplot(data_box, 'Labels', {'Lat err [m]', '|vz| [m/s]', '|att| [deg]'});
    title('Dispersion Summary'); grid on;

    sgtitle(sprintf('Campaign 2 — Touchdown Performance  |  N=%d  |  [%d,%d]m  |  wind \\sigma=%.0f m/s', ...
            sum(ok|R.status==0), h_min, h_max, sigma_vento), ...
            'FontSize', 12, 'FontWeight', 'bold');

    % ─────────────────────────────────────────────────────────────────────
    % FIGURE 2: Advanced States, Lateral Velocity and Control Effort
    % ─────────────────────────────────────────────────────────────────────
    figure('Name','Campaign 2 — Dynamics and Control Effort','NumberTitle','off',...
           'Position',[80 80 1400 700]);

    % 1. Lateral Velocity histogram (tip-over risk)
    % Velocidade lateral no toque — indicador de risco de tombamento.
    subplot(2,3,1);
    histogram(vlat, 20, 'FaceColor', [0.6 0.4 0.8]); hold on;
    xline(prctile(vlat,95), 'r--', sprintf('p95=%.2fm/s',prctile(vlat,95)), 'LineWidth', 1.5);
    xline(1.5, 'k--', 'Req (1.0)');
    xlabel('Lateral Velocity v_{lat} [m/s]'); ylabel('Count');
    title('Lateral Velocity at TD'); grid on;

    % 2. Correlation: Lateral Error vs Lateral Velocity
    % Erro de posição vs erro de velocidade lateral, cor = atitude.
    subplot(2,3,2);
    scatter(el, vlat, 35, att, 'filled'); hold on;
    colorbar; colormap(gca,'parula');
    yline(1.0, 'k--'); xline(5.0, 'k--');
    xlabel('Lateral Error [m]'); ylabel('Lateral Velocity [m/s]');
    title('Position Error vs Velocity Error (color = attitude°)'); grid on;

    % 3. Burn Time histogram (direct proxy for propellant mass)
    % Tempo de burn ≈ proxy do propelente gasto (esforço de controle).
    subplot(2,3,3);
    t_burn_ok = R.ign.T_burn(ok);
    histogram(t_burn_ok, 20, 'FaceColor', [0.4 0.7 0.7]); hold on;
    xline(median(t_burn_ok), 'b--', sprintf('Med=%.1fs',median(t_burn_ok)), 'LineWidth', 1.5);
    xlabel('Burn Time T_{burn} [s]'); ylabel('Count');
    title('Control Effort / Propellant Consumption'); grid on;

    % 4. How the S-curve reacts to disturbances: h_ign vs vz_ign
    % Decisão da S-curve: altitude vs velocidade na ignição, cor = erro lateral.
    subplot(2,3,4);
    h_ign_ok = R.ign.h(ok);
    vz_ign_ok = R.ign.vz(ok);
    scatter(vz_ign_ok, h_ign_ok, 35, el, 'filled'); hold on;
    colorbar; colormap(gca,'hot');
    xlabel('Vertical Velocity at Ignition [m/s]'); ylabel('Ignition Altitude h_{ign} [m]');
    title('S-curve Decision (color = lat\_err)'); grid on;

    % 5. Angular dynamics: theta_max during flight
    % Pico de theta durante o voo — margem de autoridade angular.
    subplot(2,3,5);
    tm_ok = R.theta_max(ok);
    histogram(tm_ok, 20, 'FaceColor', [0.8 0.6 0.2]); hold on;
    xline(prctile(tm_ok,95), 'r--', sprintf('p95=%.1f°',prctile(tm_ok,95)), 'LineWidth', 1.5);
    xlabel('Maximum \theta during flight [deg]'); ylabel('Count');
    title('Authority and Maximum Angular Deviation'); grid on;

    % 6. 3D scatter of touchdown velocities
    % Espaço de estados das velocidades finais (vx, vy, vz), cor = v_lat.
    subplot(2,3,6);
    scatter3(R.vx_td(ok), R.vy_td(ok), R.vz_td(ok), 30, vlat, 'filled');
    xlabel('v_x [m/s]'); ylabel('v_y [m/s]'); zlabel('v_z [m/s]');
    title('Final Velocity State Space');
    grid on; view(45,30); colormap(gca,'cool');

    sgtitle('Campaign 2 — Advanced Dynamics, Actuators and Margins', ...
            'FontSize', 12, 'FontWeight', 'bold');
end
% =========================================================================
% salvar_crashes — exporta todas as runs com status==0 para CSV
%
% Uso (a partir do main, depois das chamadas plotar_*):
%   salvar_crashes(R,  'crashes_envelope.csv');
%   salvar_crashes(R2, 'crashes_dense.csv');
%
% Cada linha do CSV é uma run que CRASHOU (chegou ao chão violando
% requisito de erro lateral OU velocidade vertical OU atitude).
% Inclui as 10 condições iniciais, 6 estados na ignição, 13 métricas
% de touchdown, e a mensagem de falha.
%
% =========================================================================
function salvar_crashes(R, filename)
% Monta uma table só com as runs que crasharam (status==0) e grava em CSV,
% para inspeção posterior das condições que levaram à falha.
    if nargin < 2, filename = 'crashes.csv'; end

    idx = find(R.status == 0);
    if isempty(idx)
        fprintf('  [salvar_crashes] Nenhum crash detectado. CSV nao gerado.\n');
        return;
    end

    T = table();
    T.run_id      = idx;
    T.faixa_idx   = R.faixa_idx(idx);
    T.h0_nominal  = R.h0_nominal(idx);

    % ── Condicoes iniciais ──────────────────────────────────────────────
    T.ic_h0       = R.ic.h0(idx);
    T.ic_vz0      = R.ic.vz0(idx);
    T.ic_x0       = R.ic.x0(idx);
    T.ic_y0       = R.ic.y0(idx);
    T.ic_phi0_deg = rad2deg(R.ic.phi0(idx));
    T.ic_theta0_deg = rad2deg(R.ic.theta0(idx));
    T.ic_q0_dps   = rad2deg(R.ic.q0(idx));
    T.ic_m0       = R.ic.m0(idx);
    T.ic_vw_x     = R.ic.vw_x(idx);
    T.ic_vw_y     = R.ic.vw_y(idx);

    % ── Estado na ignicao ───────────────────────────────────────────────
    T.ign_h       = R.ign.h(idx);
    T.ign_vz      = R.ign.vz(idx);
    T.ign_x       = R.ign.x(idx);
    T.ign_y       = R.ign.y(idx);
    T.ign_t       = R.ign.t_ign(idx);
    T.ign_T_burn  = R.ign.T_burn(idx);

    % ── Touchdown ──────────────────────────────────────────────────────
    T.td_x        = R.x_td(idx);
    T.td_y        = R.y_td(idx);
    T.td_z        = R.z_td(idx);
    T.td_vx       = R.vx_td(idx);
    T.td_vy       = R.vy_td(idx);
    T.td_vz       = R.vz_td(idx);
    T.td_phi_deg  = rad2deg(R.phi_td(idx));
    T.td_theta_deg= rad2deg(R.theta_td(idx));
    T.td_psi_deg  = rad2deg(R.psi_td(idx));
    T.td_err_lat  = R.err_lat(idx);
    T.td_t_sim    = R.t_sim(idx);
    T.td_theta_max_deg = R.theta_max(idx);
    T.td_phi_max_deg   = R.phi_max(idx);

    % ── Mensagem (motivo da falha) ─────────────────────────────────────
    T.msg = string(R.msg(idx));

    writetable(T, filename);

    fprintf('  [salvar_crashes] %d crash(es) salvo(s) em %s\n', ...
            numel(idx), filename);
end
% #########################################################################
%  FUNÇÕES: S-CURVE (auto-contida, sem Simulink, 1000 runs com drop test
% no simulink demoraria muito mais que as 8 horas usando o analítico)
% #########################################################################

% =========================================================================
%  gerar_scurve_v2 — S-curve completa (drop analítico + scan + burn)
% =========================================================================
function [traj, ign] = gerar_scurve_v2( ...
    h0, vz0, x0, y0, phi0, theta0, ...
    hf,vzf,azf, xf,vxf,axf, yf,vyf,ayf, ...
    empuxo_max, T_min_motor, m0, margem, dt, g, ...
    SCAN_DT, TVIOL, ALPHA)
% Gera a trajetória de referência: queda livre analítica até a ignição,
% varredura da janela de ignição viável e S-curve (quíntico) do burn.
% Lança erro se a manobra é inviável para estas CIs (capturado pelo chamador).

    ign = struct('h_ign',NaN,'vz_ign',NaN,'x_ign',NaN,'y_ign',NaN,...
                 'T_burn',NaN,'t_ign',NaN);
    psi0 = 0;

    % ── Corrigir velocidades inerciais pela atitude ──────────────────
    % Se há inclinação inicial, a velocidade vertical do corpo gera
    % componentes laterais no inercial; rotaciona vz0 por R_IB para obtê-las.
    vx0 = 0; vy0 = 0;
    if abs(phi0) > 1e-6 || abs(theta0) > 1e-6
        R_IB = rot_zyx(phi0, theta0, psi0);
        v_I  = R_IB * [0; 0; vz0];
        vx0 = v_I(1);  vy0 = v_I(2);
    end

    % ── Drop test analítico ──────────────────────────────────────────
    % Tempo de impacto da queda livre (raiz positiva da quadrática em t).
    disc = vz0^2 + 2*g*h0;
    if disc < 0, error('mc:no_impact','Sem impacto: disc<0'); end
    t_imp = (-vz0 + sqrt(disc)) / g;

    % Perfis analíticos da queda (atitude congelada, sem empuxo).
    tv = (0:0.01:t_imp)';
    drop.t  = tv;
    drop.h  = h0 + vz0*tv - 0.5*g*tv.^2;
    drop.vz = vz0 - g*tv;
    drop.x  = x0 + vx0*tv;
    drop.y  = y0 + vy0*tv;
    drop.vx = vx0*ones(size(tv));
    drop.vy = vy0*ones(size(tv));
    drop.phi   = phi0*ones(size(tv));
    drop.theta = theta0*ones(size(tv));
    drop.psi   = zeros(size(tv));

    N_drop = length(tv);

    % ── Scan ─────────────────────────────────────────────────────────
    % Subamostra a queda no passo SCAN_DT e testa cada candidato a ignição.
    dt_d = mean(diff(tv));
    skip = max(1, round(SCAN_DT / dt_d));
    sidx = 1:skip:N_drop;
    Ns   = length(sidx);

    % Resultados por candidato: tempo, altitude, h_min, violação, validade.
    st = NaN(Ns,1); sh = NaN(Ns,1); shm = NaN(Ns,1);
    stv = NaN(Ns,1); sv = false(Ns,1);

    for k = 1:Ns
        ii = sidx(k);
        t_i=drop.t(ii); h_i=drop.h(ii);
        x_i=drop.x(ii); y_i=drop.y(ii);
        vx_i=drop.vx(ii); vy_i=drop.vy(ii); vz_i=drop.vz(ii);
        phi_i=drop.phi(ii); theta_i=drop.theta(ii);

        st(k) = t_i; sh(k) = h_i;
        if h_i <= hf + 0.5, continue; end   % perto do solo: ignição inviável

        % Aceleração no instante (T_min projetado pela atitude congelada).
        ax_i= (T_min_motor/m0)*sin(theta_i)*cos(phi_i);
        ay_i=-(T_min_motor/m0)*sin(phi_i);
        az_i= (T_min_motor/m0)*cos(theta_i)*cos(phi_i) - g;

        % Testa a viabilidade da S-curve ignitando aqui.
        [ok_k, hm_k, tv_k] = scan_rapido( ...
            h_i,vz_i,az_i,hf,vzf,azf, ...
            x_i,vx_i,ax_i,xf,vxf,axf, ...
            y_i,vy_i,ay_i,yf,vyf,ayf, ...
            empuxo_max,T_min_motor,m0,margem,dt,g);
        sv(k)=ok_k; shm(k)=hm_k; stv(k)=tv_k;
    end

    % ── Janela ───────────────────────────────────────────────────────
    % Limite tardio: último candidato viável que não fura o chão.
    il = find(sv & shm >= -0.5);
    if isempty(il), error('mc:no_late','Sem frame com h_min>=0'); end
    t_late = st(il(end));

    % Limite precoce: primeiro candidato sem violar T_min (ou o de menor violação).
    ie = find(sv & stv <= TVIOL);
    if ~isempty(ie)
        t_early = st(ie(1));
    else
        [~,ib] = min(stv);
        t_early = st(ib);
    end

    % Janela invertida = não há instante de ignição viável.
    if t_late < t_early
        error('mc:env_fechado','Envelope fechado');
    end

    % ── Seleção ──────────────────────────────────────────────────────
    % Ignita em ALPHA da janela (0.5 = meio) e ajusta ao ponto de grade.
    t_ign = t_late - ALPHA*(t_late - t_early);
    [~, ii] = min(abs(tv - t_ign));
    t_ign = tv(ii);

    % Estado de queda livre exatamente no instante de ignição escolhido.
    h_ign = drop.h(ii); vz_ign = drop.vz(ii);
    x_ign = drop.x(ii); y_ign = drop.y(ii);
    vx_ign = drop.vx(ii); vy_ign = drop.vy(ii);
    phi_ign = drop.phi(ii); theta_ign = drop.theta(ii);

    ax_ign= (T_min_motor/m0)*sin(theta_ign)*cos(phi_ign);
    ay_ign=-(T_min_motor/m0)*sin(phi_ign);
    az_ign= (T_min_motor/m0)*cos(theta_ign)*cos(phi_ign)-g;

    ign.h_ign=h_ign; ign.vz_ign=vz_ign; ign.x_ign=x_ign; ign.y_ign=y_ign;
    ign.t_ign=t_ign;

    % ── Burn ─────────────────────────────────────────────────────────
    % S-curve definitiva (quíntico nos 3 canais) do estado de ignição ao toque.
    tb = calc_3d(h_ign,vz_ign,az_ign,hf,vzf,azf, ...
                 x_ign,vx_ign,ax_ign,xf,vxf,axf, ...
                 y_ign,vy_ign,ay_ign,yf,vyf,ayf, ...
                 empuxo_max,T_min_motor,m0,margem,dt,g);
    ign.T_burn = tb.T;

    % ── Concatenar ───────────────────────────────────────────────────
    % Junta o trecho de queda (até a ignição) com o trecho de burn,
    % deslocando o tempo do burn para começar em t_ign.
    if ii > 1, is = 1:(ii-1); else, is = []; end
    Nf = length(is);

    ts = tb.time + t_ign;
    traj.time    = [tv(is); ts];
    traj.h_d     = [drop.h(is); tb.h_d];
    traj.vz_d    = [drop.vz(is); tb.vz_d];
    traj.az_d    = [-g*ones(Nf,1); tb.az_d];   % na queda az = -g
    traj.x_d     = [drop.x(is); tb.x_d];
    traj.vx_d    = [drop.vx(is); tb.vx_d];
    traj.ax_d    = [zeros(Nf,1); tb.ax_d];     % na queda não há acel. lateral
    traj.y_d     = [drop.y(is); tb.y_d];
    traj.vy_d    = [drop.vy(is); tb.vy_d];
    traj.ay_d    = [zeros(Nf,1); tb.ay_d];
    traj.theta_d = [drop.theta(is); tb.theta_d];
    traj.phi_d   = [drop.phi(is); tb.phi_d];
    traj.Tff_d   = [zeros(Nf,1); tb.Tff_d];

    psi_t = [drop.psi(is); zeros(size(ts))];

    % Body velocities
    % Converte a velocidade de referência inercial → corpo (R_BI = R'),
    % ponto a ponto, pois o controlador consome velocidade no corpo.
    u_t=zeros(size(traj.time));
    v_t=zeros(size(traj.time));
    w_t=zeros(size(traj.time));
    for i=1:length(traj.time)
        R=rot_zyx(traj.phi_d(i),traj.theta_d(i),psi_t(i));
        vb=R'*[traj.vx_d(i);traj.vy_d(i);traj.vz_d(i)];
        u_t(i)=vb(1); v_t(i)=vb(2); w_t(i)=vb(3);
    end
    traj.u_d=u_t; traj.v_d=v_t; traj.w_d=w_t;

    traj.T=t_ign+tb.T; traj.T_burn=tb.T; traj.t_ign=t_ign; traj.dt=dt;

    % timeseries
    % Empacota cada perfil como timeseries para os blocos From Workspace.
    tt = traj.time;
    traj.h_ref     = timeseries(traj.h_d,tt,'Name','h_ref');
    traj.vz_ref    = timeseries(traj.vz_d,tt,'Name','vz_ref');
    traj.az_ref    = timeseries(traj.az_d,tt,'Name','az_ref');
    traj.x_ref     = timeseries(traj.x_d,tt,'Name','x_ref');
    traj.vx_ref    = timeseries(traj.vx_d,tt,'Name','vx_ref');
    traj.ax_ref    = timeseries(traj.ax_d,tt,'Name','ax_ref');
    traj.y_ref     = timeseries(traj.y_d,tt,'Name','y_ref');
    traj.vy_ref    = timeseries(traj.vy_d,tt,'Name','vy_ref');
    traj.ay_ref    = timeseries(traj.ay_d,tt,'Name','ay_ref');
    traj.theta_ref = timeseries(traj.theta_d,tt,'Name','theta_ref');
    traj.phi_ref   = timeseries(traj.phi_d,tt,'Name','phi_ref');
    traj.Tff_ref   = timeseries(traj.Tff_d,tt,'Name','Tff_ref');
    traj.u_ref     = timeseries(u_t,tt,'Name','u_ref');
    traj.v_ref     = timeseries(v_t,tt,'Name','v_ref');
    traj.w_ref     = timeseries(w_t,tt,'Name','w_ref');
end

% =========================================================================
%  calc_3d — S-curve acoplada 3D (burn)
% =========================================================================
function tb = calc_3d(h0,vz0,az0,hf,vzf,azf,...
    x0,vx0,ax0,xf,vxf,axf,...
    y0,vy0,ay0,yf,vyf,ayf,...
    empuxo_max,~,m0,margem,dt,g)
% Resolve a S-curve do burn: dimensiona o tempo T para o pico de empuxo
% bater na margem do empuxo máximo (fzero) e devolve os perfis dos 3 canais.

    dh = h0-hf;
    % Função-resíduo: pico de empuxo(T) - empuxo permitido = 0.
    f=@(T_) peak3d(T_,h0,vz0,az0,hf,vzf,azf,dh,dt,...
        x0,vx0,ax0,xf,vxf,axf,y0,vy0,ay0,yf,vyf,ayf,m0,g)-(margem*empuxo_max);
    T_ot = fzero(f,[5,300],optimset('TolX',1e-4,'Display','off'));

    % Perfis do canal vertical (rz) e dos dois laterais (rc) para o T achado.
    [~,t,~,~,~,~,h_d,vz_d,az_d]=rz(T_ot,h0,vz0,az0,hf,vzf,azf,dh,dt);
    [~,x_d,vx_d,ax_d]=rc(T_ot,x0,vx0,ax0,xf,vxf,axf,t);
    [~,y_d,vy_d,ay_d]=rc(T_ot,y0,vy0,ay0,yf,vyf,ayf,t);

    tb.T=T_ot; tb.time=t; tb.dt=dt;
    tb.h_d=h_d; tb.vz_d=vz_d; tb.az_d=az_d;
    tb.x_d=x_d; tb.vx_d=vx_d; tb.ax_d=ax_d;
    tb.y_d=y_d; tb.vy_d=vy_d; tb.ay_d=ay_d;
    % Atitude de referência: inclinação que gera as acelerações laterais.
    tb.theta_d = atan2(ax_d, g+az_d);
    tb.phi_d   = -atan2(ay_d, g+az_d);
    tb.Tff_d   = m0*sqrt(ax_d.^2+ay_d.^2+(g+az_d).^2);   % empuxo feedforward
end

% =========================================================================
%  scan_rapido
% =========================================================================
function [ok,hm,tv] = scan_rapido(h0,vz0,az0,hf,vzf,azf,...
    x0,vx0,ax0,xf,vxf,axf,y0,vy0,ay0,yf,vyf,ayf,...
    empuxo_max,T_min_motor,m0,margem,dt,g)
% Versão leve de calc_3d usada no scan: dimensiona T, mede a altitude mínima
% (hm) e a fração de tempo (tv) em que o empuxo cairia abaixo de T_min.
% Devolve ok=false silenciosamente se o solver falhar (candidato descartado).
    ok=false; hm=NaN; tv=NaN;
    dh=h0-hf; if dh<=0,return;end
    try
        f=@(T_)peak3d(T_,h0,vz0,az0,hf,vzf,azf,dh,dt,...
            x0,vx0,ax0,xf,vxf,axf,y0,vy0,ay0,yf,vyf,ayf,...
            m0,g)-(margem*empuxo_max);
        T_ot=fzero(f,[5,300],optimset('TolX',1e-3,'Display','off'));
    catch, return; end
    try [~,t,~,~,~,~,h_d,~,az_d]=rz(T_ot,h0,vz0,az0,hf,vzf,azf,dh,dt);
    catch, return; end
    hm=min(h_d); ok=true;
    if hm<0, tv=1.0; return; end   % furou o chão: viola tudo
    try
        [~,~,~,ax_d]=rc(T_ot,x0,vx0,ax0,xf,vxf,axf,t);
        [~,~,~,ay_d]=rc(T_ot,y0,vy0,ay0,yf,vyf,ayf,t);
    catch, tv=1.0; return; end
    Tf=m0*sqrt(ax_d.^2+ay_d.^2+(g+az_d).^2);
    tv=sum(Tf<T_min_motor)/length(Tf);
end

% =========================================================================
%  peak3d — pico de empuxo 3D
% =========================================================================
function Tp=peak3d(T_,h0,vz0,az0,hf,vzf,azf,dh,dt,...
    x0,vx0,ax0,xf,vxf,axf,y0,vy0,ay0,yf,vyf,ayf,m0,g)
% Pico de empuxo feedforward para um tempo de burn T: resolve os 3 canais e
% retorna max(Tff). Usada como função-resíduo do fzero em calc_3d/scan_rapido.
    try
        [~,t,~,~,~,~,~,~,az_d]=rz(T_,h0,vz0,az0,hf,vzf,azf,dh,dt);
        [~,~,~,ax_d]=rc(T_,x0,vx0,ax0,xf,vxf,axf,t);
        [~,~,~,ay_d]=rc(T_,y0,vy0,ay0,yf,vyf,ayf,t);
        Tp=max(m0*sqrt(ax_d.^2+ay_d.^2+(g+az_d).^2));
    catch, Tp=NaN; end
end

% =========================================================================
%  rz — resolver S-curve eixo Z
% =========================================================================
function [c,t,tau,s,sd,sdd,h,v,a]=rz(T_,h0,vz0,az0,hf,vzf,azf,dh,dt)
% Polinômio quíntico do canal vertical em variável normalizada tau=t/T.
% M é a matriz de contorno (valor/1a/2a derivada em tau=0 e 1); b traz as
% condições já normalizadas por T e dh. Retorna s e desnormaliza para h/v/a.
    M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 2 0 0 0;
       1 1 1 1 1 1;0 1 2 3 4 5;0 0 2 6 12 20];
    b=[0;-vz0*T_/dh;-az0*T_^2/dh;1;-vzf*T_/dh;-azf*T_^2/dh];
    c=M\b; t=(0:dt:T_)'; tau=t/T_;
    s=c(1)+c(2)*tau+c(3)*tau.^2+c(4)*tau.^3+c(5)*tau.^4+c(6)*tau.^5;
    sd=c(2)+2*c(3)*tau+3*c(4)*tau.^2+4*c(5)*tau.^3+5*c(6)*tau.^4;
    sdd=2*c(3)+6*c(4)*tau+12*c(5)*tau.^2+20*c(6)*tau.^3;
    h=h0-dh*s; v=-(dh/T_)*sd; a=-(dh/T_^2)*sdd;
end

% =========================================================================
%  rc — resolver canal lateral
% =========================================================================
function [c,p,v,a]=rc(T_,p0,v0,a0,pf,vf,af,t)
% Polinômio quíntico de um canal lateral (x ou y), direto em unidades
% físicas. v0/a0 entram escalados por T e T^2 porque as derivadas são em tau.
    tau=t/T_;
    M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 2 0 0 0;
       1 1 1 1 1 1;0 1 2 3 4 5;0 0 2 6 12 20];
    b=[p0;T_*v0;T_^2*a0;pf;T_*vf;T_^2*af];
    c=M\b;
    p=c(1)+c(2)*tau+c(3)*tau.^2+c(4)*tau.^3+c(5)*tau.^4+c(6)*tau.^5;
    v=(c(2)+2*c(3)*tau+3*c(4)*tau.^2+4*c(5)*tau.^3+5*c(6)*tau.^4)/T_;
    a=(2*c(3)+6*c(4)*tau+12*c(5)*tau.^2+20*c(6)*tau.^3)/T_^2;
end

% =========================================================================
%  rot_zyx — rotação ZYX
% =========================================================================
function R=rot_zyx(phi,theta,psi)
% Matriz de rotação corpo→inercial (sequência ZYX), versão local equivalente
% a rotation_matrix.m, para manter a S-curve auto-contida.
    cp=cos(phi);sp=sin(phi);ct=cos(theta);st=sin(theta);
    cs=cos(psi);ss=sin(psi);
    R=[ct*cs,sp*st*cs-cp*ss,cp*st*cs+sp*ss;
       ct*ss,sp*st*ss+cp*cs,cp*st*ss-sp*cs;
       -st,sp*ct,cp*ct];
end

% =========================================================================
%  Utilitários
% =========================================================================
% tg: ruído gaussiano truncado em ±3σ (rejeita amostras na cauda extrema).
function v=tg(s), v=s*randn(); while abs(v)>3*s, v=s*randn(); end, end

% chk: imprime uma linha da tabela de requisitos (OK se valor <= limite).
function chk(n,v,r)
    if v<=r, s='  OK'; else, s='  FALHA'; end
    fprintf('║  %-22s  %7.2f  %7.2f  %s\n',n,v,r,s);
end

% iff: seleção condicional inline (operador ternário inexistente em MATLAB).
function r=iff(c,a,b), if c,r=a;else,r=b;end, end

% print_progress_v2: imprime a linha de progresso de uma run, com ETA.
function print_progress_v2(k,N,h0,status,el,vz,vlat,tm,dt_r,t0)
    te=toc(t0); ta=te/k; eta=ta*(N-k); pct=100*k/N;
    tags={'SIM-FL','S-FAIL',' CRASH','  ???  ','  OK  '};
    tag=tags{min(5,max(1,status+3))};   % mapeia status (-2..1) para rótulo
    if isnan(el)
        fprintf('  %3d  | %5.0f  |  %s  |    ---    |   ---   |   ---   |    ---    | %5.1fs | %s (%.0f%%)\n',...
                k,h0,tag,dt_r,feta(eta),pct);
    else
        fprintf('  %3d  | %5.0f  |  %s  | %7.2f m | %+5.1f   | %5.2f   | %6.1f°   | %5.1fs | %s (%.0f%%)\n',...
                k,h0,tag,el,vz,vlat,tm,dt_r,feta(eta),pct);
    end
end

% feta: formata segundos como "~Ns" ou "~Nmin" para o display do ETA.
function s=feta(sec)
    if sec<60, s=sprintf('~%.0fs',sec);
    else, s=sprintf('~%.0fmin',sec/60); end
end