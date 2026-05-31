% =========================================================================
%  RUN_SCURVE.M
%  Geracao de trajetoria de referencia para pouso propulsivo (suicide burn)
%  de um foguete reutilizavel, em tres eixos (3D), com selecao automatica
%  do instante de ignicao do motor.
% -------------------------------------------------------------------------
%  Geracao do perfil de referencia
%  (posicao, velocidade, aceleracao e atitude desejadas ao longo do tempo)
%  que sera posteriormente seguido pela malha de controle no modelo
%  dinamico do foguete. Esse perfil e o que comumente se chama de
%  "guidance" ou trajetoria nominal de feedforward.

%  Solucao adotada ("Predictive True-State Scan"):
%
%      1. DROP TEST - Simula a queda livre no proprio modelo Simulink
%         (planta_foguete.slx) com o motor desligado (empuxo e TVC nulos)
%         ate o foguete atingir h ~ 0. Dessa simulacao extraem-se as
%         series temporais REAIS de h, x, y, vx, vy, vz, phi e theta.
%
%      2. TRUE-STATE SCAN - Percorre o vetor de tempo da queda (a cada
%         ~0.2 s). Em cada instante candidato usa o estado REAL como CI
%         da curva de pouso e registra:
%           h_min  - altitude minima do polinomio gerado (verifica se a
%                    trajetoria fura o solo, o chamado "toupeira check");
%           T_viol - fracao do burn em que o empuxo exigido fica abaixo
%                    do empuxo minimo do motor.
%
%      3. JANELA REAL - A partir do scan determina:
%           t_late  - maior instante (mais proximo do solo) em que ainda
%                     ha frenagem viavel sem furar o solo (h_min >= 0);
%           t_early - menor instante (mais alto) em que a violacao de
%                     empuxo minimo fica dentro da tolerancia.
%
%      4. SELECAO POR alpha_ign - O usuario escolhe um valor entre 0 e 1,
%         de agressivo (alpha = 0, ignicao tardia em t_late) a conservador
%         (alpha = 1, ignicao precoce em t_early). O script entao extrai o
%         estado exato naquele instante e gera a trajetoria final.
% -------------------------------------------------------------------------
%  REQUISITOS DO MODELO SIMULINK
%    O modelo 'planta_foguete.slx' precisa estar no path do MATLAB.
%  SAIDAS PRINCIPAIS (geradas no workspace ao final da execucao)
%    traj     - struct com a trajetoria de referencia completa (queda +
%               burn), incluindo timeseries prontas para o Simulink.
%    envelope - struct com os dados do drop test e do scan, util para
%               diagnostico e para os graficos.
%    drop     - struct com as series temporais brutas do drop test.
% =========================================================================

% Carrega o conjunto de parametros fisicos e operacionais do foguete
% (massa, empuxo maximo, gravidade, passo de integracao, condicoes de
% contorno, etc.). Toda a configuracao numerica vem desta struct 'p'.
p = falcon9_params();

% ── Condicoes de contorno verticais (eixo Z / altitude) ──────────────────
% Indice 0 = estado inicial da trajetoria; indice f = estado final (pouso).
% h: altitude, vz: velocidade vertical, az: aceleracao vertical.
h0  = p.h0;    vz0 = p.vz0;   az0 = p.az0;
hf  = p.hf;    vzf = p.vzf;   azf = p.azf;

% ── Condicoes de contorno laterais no eixo X ─────────────────────────────
% Usa local_get para ler o campo se existir ou assumir 0 como padrao,
% tornando o script tolerante a structs de parametros incompletas.
x0  = local_get(p, 'x0',  0);  vx0 = local_get(p, 'vx0', 0);
ax0 = local_get(p, 'ax0', 0);  xf  = local_get(p, 'xf',  0);
vxf = local_get(p, 'vxf', 0);  axf = local_get(p, 'axf', 0);

% ── Condicoes de contorno laterais no eixo Y ─────────────────────────────
y0  = local_get(p, 'y0',  0);  vy0 = local_get(p, 'vy0', 0);
ay0 = local_get(p, 'ay0', 0);  yf  = local_get(p, 'yf',  0);
vyf = local_get(p, 'vyf', 0);  ayf = local_get(p, 'ayf', 0);

% ── Atitude inicial do foguete (angulos de Euler) ────────────────────────
phi_0   = local_get(p, 'phi0',   0);
theta_0 = local_get(p, 'theta0', 0);
psi_0   = local_get(p, 'psi0',   0);

% ── Projecao das velocidades do corpo para o referencial inercial ────────
% As velocidades p.vx0/vy0/vz0 sao fornecidas no referencial do CORPO
% (consistente com a inicializacao da simulacao). O polinomio quintico,
% no entanto, opera em coordenadas inerciais. Por isso, quando ha atitude
% inicial nao nula, projetamos a CI de velocidade para o inercial antes
% de usa-la. Se a atitude inicial e essencialmente nula, a projecao seria
% identidade e o bloco e simplesmente ignorado.
if abs(phi_0) > 1e-6 || abs(theta_0) > 1e-6 || abs(psi_0) > 1e-6
    R_IB = rotation_matrix(phi_0, theta_0, psi_0);
    v_I  = R_IB * [vx0; vy0; vz0];      % usa as 3 componentes de velocidade
    vx0  = v_I(1);
    vy0  = v_I(2);
    vz0  = v_I(3);                       % atualiza tambem vz0 (era w no corpo)
    fprintf('  Projecao body->inercial: (vx,vy,vz) = (%.3f, %.3f, %.3f) m/s\n', ...
            vx0, vy0, vz0);
end

% ── Parametros fisicos e operacionais usados na otimizacao ───────────────
empuxo_max   = p.T_max;          % empuxo maximo do motor [N]
throttle_min = p.throttle_min;   % fracao minima de acionamento do motor
m0_pouso     = p.m0;             % massa do foguete na fase de pouso [kg]
margem       = 0.90;             % margem de seguranca sobre o empuxo maximo
dt           = p.dt;             % passo de tempo da discretizacao [s]
g            = p.g;              % aceleracao da gravidade [m/s^2]

% Empuxo minimo efetivo: o motor nao consegue produzir menos do que isso
% quando ligado (limite inferior de garganta / throttle).
T_min_motor  = empuxo_max * throttle_min;

% ── Condicoes finais de pouso ────────────────────────────────────────────
% (Mantidas como definidas acima; consumidas pelo scurve_calc_3d.)

% ── CONFIGURACOES DEFINIDAS PELO USUARIO ─────────────────────────────────
alpha_ign_user = 0.5;    % 0 = agressivo (ignicao tardia); 1 = conservador (ignicao precoce)
SLX_MODEL      = 'planta_foguete';    % nome do modelo Simulink a ser simulado
T_SIM_LIMITE   = 15;                  % tempo limite do drop test [s]
SCAN_DT_RES    = 0.2;                 % resolucao temporal desejada para o scan [s]
T_VIOL_TOL     = 0.03;                % tolerancia de violacao de empuxo minimo (3%)

fprintf('\n=========================================================\n');
fprintf('  RUN_SCURVE 3D — v6: Predictive True-State Scan\n');
fprintf('=========================================================\n');

% =========================================================================
%  DIRETRIZ 1 — DROP TEST: simulacao preditiva da queda livre
% -------------------------------------------------------------------------
%  Roda o modelo planta_foguete.slx com o motor desligado (empuxo e TVC
%  nulos) e captura o estado real do foguete ao longo do tempo, incluindo
%  os efeitos aerodinamicos e o desvio inercial decorrente da atitude.
% =========================================================================

% ── Sinaliza ao Simulink para cortar o motor durante a queda ─────────────
% Este e o unico assignin necessario; o modelo le esta flag e zera o
% empuxo e o controle de vetorizacao de empuxo (TVC).
assignin('base', 'FLAG_MOTOR_OFF', true);

% ── Executa a simulacao da queda (sem forcar condicoes iniciais) ─────────
fprintf('    Simulando...');
t_sim_start = tic;

try
    sim_out = sim(SLX_MODEL, [0, T_SIM_LIMITE]);
catch ME
    % Caso a simulacao falhe, interrompe com mensagem clara em vez de
    % deixar o erro propagar de forma confusa.
    error('run_scurve:drop_test_fail', 'Drop test falhou: %s', ME.message);
end

t_sim_elapsed = toc(t_sim_start);
fprintf(' concluído em %.2f s\n', t_sim_elapsed);

% ── Extracao robusta dos dados logados ───────────────────────────────────
% O procedimento abaixo e identico ao usado em metrica_qualidade_pouso.m,
% garantindo consistencia entre os scripts. O sinal logado vem como uma
% timeseries cujo campo Data pode ter dimensoes diferentes dependendo da
% configuracao do modelo; por isso normalizamos o formato para uma matriz
% [amostras x estados].
ts = sim_out.get('x_sim_planta_pura');
drop.t = ts.Time;

dados_raw = ts.Data;
if ndims(dados_raw) == 3
    % Formato 3D tipico de logsout: remove dimensoes unitarias e transpoe.
    dados = squeeze(dados_raw)';
else
    if size(dados_raw, 1) == length(drop.t)
        % Ja esta no formato [amostras x estados].
        dados = dados_raw;
    else
        % Esta transposto; corrige.
        dados = dados_raw';
    end
end

% ── Mapeamento das colunas de estado e correcao de referencial ───────────
% Colunas 1-3: posicao inercial (x, y, h).
drop.x     = dados(:, 1);
drop.y     = dados(:, 2);
drop.h     = dados(:, 3);

% Colunas 4-6: velocidades expressas no REFERENCIAL DO CORPO (u, v, w).
% Atencao: o Simulink devolve velocidades no corpo, nao no inercial.
u_B = dados(:, 4);
v_B = dados(:, 5);
w_B = dados(:, 6);

% Colunas 7-9: angulos de atitude
drop.phi   = dados(:, 7);
drop.theta = dados(:, 8);
drop.psi   = dados(:, 9);

% Prealoca os vetores de velocidade ja no referencial inercial.
drop.vx = zeros(size(u_B));
drop.vy = zeros(size(v_B));
drop.vz = zeros(size(w_B));

% ── Rotaciona as velocidades do corpo para o inercial, quadro a quadro ───
% Como a atitude varia ao longo da queda, a matriz de rotacao precisa ser
% recalculada em cada instante de tempo.
for i = 1:length(drop.t)
    phi_i = drop.phi(i);
    th_i  = drop.theta(i);
    psi_i = drop.psi(i);

    % Matriz de rotacao R_IB (corpo -> inercial)
    % (sequencia de Euler 3-2-1).
    R_IB = [cos(th_i)*cos(psi_i), sin(phi_i)*sin(th_i)*cos(psi_i) - cos(phi_i)*sin(psi_i), cos(phi_i)*sin(th_i)*cos(psi_i) + sin(phi_i)*sin(psi_i);
            cos(th_i)*sin(psi_i), sin(phi_i)*sin(th_i)*sin(psi_i) + cos(phi_i)*cos(psi_i), cos(phi_i)*sin(th_i)*sin(psi_i) - sin(phi_i)*cos(psi_i);
           -sin(th_i),            sin(phi_i)*cos(th_i),                                  cos(phi_i)*cos(th_i)];

    V_I = R_IB * [u_B(i); v_B(i); w_B(i)];

    drop.vx(i) = V_I(1);
    drop.vy(i) = V_I(2);
    drop.vz(i) = V_I(3);
end

% ── Trunca as series no instante do impacto (h ~ 0) ──────────────────────
% Procura o primeiro indice em que a altitude cruza o solo.
idx_solo = find(drop.h <= 0, 1, 'first');
if isempty(idx_solo)
    % Se o foguete nao chegou ao solo dentro do tempo limite, usa toda a
    % serie disponivel e avisa o usuario.
    idx_solo = length(drop.t);
    fprintf('    [INFO] Foguete não atingiu o chão em %.1f s (h_final=%.1f m)\n', ...
            T_SIM_LIMITE, drop.h(end));
end

% Aplica o truncamento a todos os campos da struct 'drop' de forma
% generica, percorrendo seus nomes de campo.
fn = fieldnames(drop);
for i = 1:numel(fn)
    drop.(fn{i}) = drop.(fn{i})(1:idx_solo);
end

N_drop = length(drop.t);     % numero de amostras uteis da queda
t_impacto = drop.t(end);     % instante aproximado do impacto

fprintf('    Drop Test: %d amostras, t_impacto ≈ %.2f s\n', N_drop, t_impacto);
fprintf('    Estado no impacto: h=%.1f m  vz=%.1f m/s\n', drop.h(end), drop.vz(end));
% Compara o desvio lateral real com o que a aproximacao analitica preveria,
% evidenciando o efeito do arrasto e da atitude.
fprintf('    Drift: dx=%.2f m  dy=%.2f m  (vs analítico: dx=%.2f  dy=%.2f)\n', ...
        drop.x(end) - x0, drop.y(end) - y0, ...
        vx0 * t_impacto, vy0 * t_impacto);

% ── Limpa a flag do workspace (motor volta a estar habilitado) ───────────
assignin('base', 'FLAG_MOTOR_OFF', false);

% =========================================================================
%  DIRETRIZ 2 — TRUE-STATE SCAN: varredura sobre o estado real
% -------------------------------------------------------------------------
%  Em vez de varrer um parametro abstrato, percorremos o vetor de tempo
%  real da queda. Para cada instante candidato t_i extraimos o estado REAL
%  do drop test e tentamos gerar um polinomio de pouso a partir dele.
% =========================================================================

fprintf('\n  DIRETRIZ 2 — True-State Scan\n');

% ── Define o passo de amostragem do scan ─────────────────────────────────
% Como o passo do drop test pode ser fino, pulamos amostras o suficiente
% para aproximar a resolucao desejada (SCAN_DT_RES) sem varrer pontos em
% excesso.
dt_drop     = mean(diff(drop.t));                     % passo medio do drop test
skip_N      = max(1, round(SCAN_DT_RES / dt_drop));   % quantas amostras pular
scan_idx    = 1 : skip_N : N_drop;                    % indices efetivamente varridos
N_scan      = length(scan_idx);

fprintf('    Resolução: %.3f s (skip=%d amostras, dt_drop=%.4f s)\n', ...
        skip_N * dt_drop, skip_N, dt_drop);
fprintf('    Pontos a varrer: %d\n', N_scan);

% ── Vetores que armazenam o resultado do scan ────────────────────────────
scan_t       = NaN(N_scan, 1);     % instante de cada frame avaliado
scan_h       = NaN(N_scan, 1);     % altitude no frame
scan_h_min   = NaN(N_scan, 1);     % altitude minima do polinomio gerado
scan_T_viol  = NaN(N_scan, 1);     % fracao do burn com empuxo abaixo do minimo
scan_valid   = false(N_scan, 1);   % indica se o frame gerou solucao valida

fprintf('    Varrendo...');
t_scan_start = tic;

for k = 1:N_scan
    ii = scan_idx(k);

    % ── Extrai o estado REAL neste frame da queda ───────────────────────
    t_i     = drop.t(ii);
    h_i     = drop.h(ii);
    x_i     = drop.x(ii);
    y_i     = drop.y(ii);
    vx_i    = drop.vx(ii);
    vy_i    = drop.vy(ii);
    vz_i    = drop.vz(ii);
    phi_i   = drop.phi(ii);
    theta_i = drop.theta(ii);

    scan_t(k) = t_i;
    scan_h(k) = h_i;

    % ── Sanidade: se ja estamos praticamente no solo, ignora o frame ────
    if h_i <= hf + 0.5
        continue;
    end

    % ── Aceleracao no instante de ignicao ───────────────────────────────
    % Decompoe o empuxo minimo do motor segundo a atitude real do frame
    % (continuidade de atitude, agora com atitude real
    % vinda do drop test). Isso fornece az0/ax0/ay0 coerentes para a CI
    % da curva de pouso.
    ax_i = (T_min_motor / m0_pouso) * sin(theta_i) * cos(phi_i);
    ay_i = -(T_min_motor / m0_pouso) * sin(phi_i);
    az_i = (T_min_motor / m0_pouso) * cos(theta_i) * cos(phi_i) - g;

    % ── Tenta gerar a curva de pouso usando este estado como CI ─────────
    [ok, h_min_k, T_viol_k] = scan_frame_rapido( ...
        h_i,  vz_i, az_i,  hf, vzf, azf, ...
        x_i,  vx_i, ax_i,  xf, vxf, axf, ...
        y_i,  vy_i, ay_i,  yf, vyf, ayf, ...
        empuxo_max, T_min_motor, m0_pouso, margem, dt, g);

    scan_valid(k)   = ok;
    scan_h_min(k)   = h_min_k;
    scan_T_viol(k)  = T_viol_k;
end

t_scan_elapsed = toc(t_scan_start);
fprintf(' concluído em %.2f s  (%d frames, %d válidos)\n', ...
        t_scan_elapsed, N_scan, sum(scan_valid));

% =========================================================================
%  DIRETRIZ 3 — DETERMINACAO DA JANELA DE IGNICAO REAL
% -------------------------------------------------------------------------
%  A partir dos resultados do scan, identifica os dois limites da janela
%  de ignicao viavel: o limite tardio (agressivo) e o limite precoce
%  (conservador).
% =========================================================================

fprintf('\n  DIRETRIZ 3 — Janela de Ignição Real\n');

h_floor = -0.5;   % o polinomio nao pode cruzar o solo; pequena folga negativa

% ── t_late (limite tardio / agressivo) ───────────────────────────────────
% Maior instante (mais proximo do solo) em que o polinomio ainda consegue
% frear sem furar o solo, isto e, com altitude minima acima de h_floor.
idx_late = find(scan_valid & scan_h_min >= h_floor);
if ~isempty(idx_late)
    i_late  = idx_late(end);          % ultimo frame seguro
    t_late  = scan_t(i_late);
    h_late  = scan_h(i_late);
else
    % Se nenhum frame e capaz de frear sem furar o solo, a manobra e
    % inviavel com estas condicoes iniciais.
    error('run_scurve:no_late_window', ...
        ['[CRITICO] Nenhum frame do drop test produz h_min >= 0.\n' ...
         '          A manobra é impossível com estas CI. Abortando.']);
end

% ── t_early (limite precoce / conservador) ───────────────────────────────
% Menor instante (mais alto) em que a violacao de empuxo minimo permanece
% dentro da tolerancia T_VIOL_TOL.
idx_early = find(scan_valid & scan_T_viol <= T_VIOL_TOL);
if ~isempty(idx_early)
    i_early = idx_early(1);           % primeiro frame aceitavel
    t_early = scan_t(i_early);
    h_early = scan_h(i_early);
else
    % Se nenhum frame respeita a tolerancia, escolhe o de menor violacao
    % e avisa o usuario.
    [~, ib]   = min(scan_T_viol);
    i_early   = ib;
    t_early   = scan_t(i_early);
    h_early   = scan_h(i_early);
    fprintf('    [AVISO] Nenhum frame com T_viol <= %.1f%%. Usando melhor: t=%.2f s (%.1f%%)\n', ...
            T_VIOL_TOL*100, t_early, scan_T_viol(i_early)*100);
end

% ── Verifica se a janela esta aberta ─────────────────────────────────────
% Para a manobra ser viavel, o limite tardio precisa ocorrer depois do
% limite precoce. Caso contrario, nao ha nenhum instante valido.
if t_late < t_early
    error('run_scurve:envelope_fechado', ...
        ['[CRITICO] Envelope FECHADO! t_late (%.2f s) < t_early (%.2f s).\n' ...
         '          A manobra é impossível com estas CI.\n' ...
         '          Sugestões: reduzir h0, aumentar T_max, ou mudar alvo.'], ...
         t_late, t_early);
end

% ── Relatorio resumido da janela encontrada ──────────────────────────────
fprintf('\n  ╔═══════════════════════════════════════════════════════════╗\n');
fprintf(  '  ║       JANELA DE IGNIÇÃO REAL (True-State)                ║\n');
fprintf(  '  ╠═══════════════════════════════════════════════════════════╣\n');
fprintf(  '  ║  t_early (conservador): %6.2f s   h = %6.1f m           ║\n', t_early, h_early);
fprintf(  '  ║  t_late  (agressivo):   %6.2f s   h = %6.1f m           ║\n', t_late,  h_late);
fprintf(  '  ║  Janela:   Δt = %.2f s   Δh = %.1f m                    ║\n', ...
                          t_late - t_early, h_early - h_late);
fprintf(  '  ║  Tolerância T_min:  %.1f%%                                ║\n', T_VIOL_TOL*100);
fprintf(  '  ║  Scan: %d frames em %.2f s                              ║\n', N_scan, t_scan_elapsed);
fprintf(  '  ╚═══════════════════════════════════════════════════════════╝\n');

% =========================================================================
%  DIRETRIZ 4 — SELECAO VIA alpha_ign E GERACAO FINAL
% -------------------------------------------------------------------------
%  Interpola entre os limites da janela conforme a preferencia do usuario:
%    alpha_ign = 0  -> agressivo  (t_late,  ignicao tardia, perto do solo)
%    alpha_ign = 1  -> conservador (t_early, ignicao precoce, mais alto)
%  A interpolacao linear e: t_ign = t_late - alpha_ign * (t_late - t_early)
% =========================================================================

fprintf('\n  DIRETRIZ 4 — Seleção e Geração Final\n');

alpha_ign = max(0, min(1, alpha_ign_user));     % limita alpha ao intervalo [0, 1]
t_ign     = t_late - alpha_ign * (t_late - t_early);

% ── Ajusta t_ign ao frame mais proximo do drop test ──────────────────────
% Como o estado real so existe nos instantes amostrados, "encaixamos"
% t_ign no instante amostrado mais proximo (snap).
[~, idx_ign_drop] = min(abs(drop.t - t_ign));
t_ign = drop.t(idx_ign_drop);     % instante exato do frame escolhido

% ── Extrai o estado REAL no instante de ignicao escolhido ────────────────
h_ign      = drop.h(idx_ign_drop);
x_ign      = drop.x(idx_ign_drop);
y_ign      = drop.y(idx_ign_drop);
vx_ign     = drop.vx(idx_ign_drop);
vy_ign     = drop.vy(idx_ign_drop);
vz_ign     = drop.vz(idx_ign_drop);
phi_ign    = drop.phi(idx_ign_drop);
theta_ign  = drop.theta(idx_ign_drop);

% ── Aceleracao no instante de ignicao (continuidade de atitude) ──────────
% Mesma decomposicao usada no scan, agora aplicada ao frame definitivo.
ax_ign = (T_min_motor / m0_pouso) * sin(theta_ign) * cos(phi_ign);
ay_ign = -(T_min_motor / m0_pouso) * sin(phi_ign);
az_ign = (T_min_motor / m0_pouso) * cos(theta_ign) * cos(phi_ign) - g;

fprintf('    alpha_ign = %.2f  →  t_ign = %.3f s  (snap ao frame %d)\n', ...
        alpha_ign, t_ign, idx_ign_drop);
fprintf('    Estado REAL na ignição:\n');
fprintf('      h  = %8.2f m     vz = %8.3f m/s\n', h_ign, vz_ign);
fprintf('      x  = %8.2f m     vx = %8.3f m/s\n', x_ign, vx_ign);
fprintf('      y  = %8.2f m     vy = %8.3f m/s\n', y_ign, vy_ign);
fprintf('      phi = %+6.2f°    theta = %+6.2f°\n', ...
        rad2deg(phi_ign), rad2deg(theta_ign));
fprintf('    Aceleração na ignição (T_min @ atitude real):\n');
fprintf('      ax = %+.4f   ay = %+.4f   az = %+.4f  m/s²\n', ax_ign, ay_ign, az_ign);

% ── Compara o estado real com a predicao analitica (mede o drift) ────────
% Util para quantificar quanto a aproximacao analitica erraria neste ponto.
h_analitico  = h0 + vz0*t_ign - 0.5*g*t_ign^2;
vz_analitico = vz0 - g*t_ign;
x_analitico  = x0 + vx0*t_ign;
y_analitico  = y0 + vy0*t_ign;

fprintf('    Drift (real - analítico):\n');
fprintf('      Δh  = %+.2f m   Δvz = %+.3f m/s\n', ...
        h_ign - h_analitico, vz_ign - vz_analitico);
fprintf('      Δx  = %+.2f m   Δy  = %+.2f m\n', ...
        x_ign - x_analitico, y_ign - y_analitico);

% =========================================================================
%  FASE 2 — CURVA DE POUSO DO BURN FINAL (otimizacao acoplada 3D)
% -------------------------------------------------------------------------
%  Gera o polinomio de referencia do periodo de motor ligado, partindo do
%  estado real de ignicao e terminando nas condicoes de pouso desejadas.
%  O tempo de queima T e otimizado para respeitar o empuxo disponivel.
% =========================================================================

fprintf('\n  Gerando S-curve final...\n');

traj_burn = scurve_calc_3d( ...
    h_ign, vz_ign, az_ign,  hf, vzf, azf, ...
    x_ign, vx_ign, ax_ign,  xf, vxf, axf, ...
    y_ign, vy_ign, ay_ign,  yf, vyf, ayf, ...
    empuxo_max, T_min_motor, m0_pouso, margem, dt, g);

T_burn = traj_burn.T;

fprintf('    T_burn (otimizado 3D): %.2f s\n', T_burn);
fprintf('    T_total: %.2f s (%.2f queda + %.2f burn)\n', ...
        t_ign + T_burn, t_ign, T_burn);

% =========================================================================
%  FASE 1 — REFERENCIA DA QUEDA LIVRE (extraida do drop test real)
% -------------------------------------------------------------------------
%  Em vez da aproximacao analitica, usamos diretamente os dados reais do
%  drop test ate o instante de ignicao. Assim a referencia ja incorpora
%  o desvio lateral, o arrasto e a atitude real durante a queda.
% =========================================================================

t_ff     = t_ign;            % renomeia para uso nos graficos (free-fall)
i_ff_end = idx_ign_drop;     % indice final da fase de queda

t_ff_vec = drop.t(1:i_ff_end);
N_ff     = length(t_ff_vec);

% Estado vertical durante a queda.
h_ff     = drop.h(1:i_ff_end);
vz_ff    = drop.vz(1:i_ff_end);
az_ff    = -g * ones(N_ff, 1);      

% Estado lateral X durante a queda.
x_ff     = drop.x(1:i_ff_end);
vx_ff    = drop.vx(1:i_ff_end);
ax_ff    = zeros(N_ff, 1);

% Estado lateral Y durante a queda.
y_ff     = drop.y(1:i_ff_end);
vy_ff    = drop.vy(1:i_ff_end);
ay_ff    = zeros(N_ff, 1);

% Atitude REAL durante a queda (nao zerada) e empuxo nulo (motor desligado).
theta_ff = drop.theta(1:i_ff_end);
phi_ff   = drop.phi(1:i_ff_end);
Tff_ff   = zeros(N_ff, 1);

% =========================================================================
%  CONCATENACAO: queda livre (drop test) + burn (curva de pouso)
% -------------------------------------------------------------------------
%  Junta as duas fases em vetores unicos, deslocando o tempo do burn para
%  comecar em t_ff. O ultimo ponto da queda e removido para nao duplicar
%  com o primeiro ponto do burn.
% =========================================================================

if N_ff > 1
    i_ff_sel = 1:(N_ff-1);           % remove o ultimo ponto (duplicado)
else
    i_ff_sel = [];
end

t_burn_shifted = traj_burn.time + t_ff;   % desloca o tempo do burn

% Vertical (Z / altitude).
t_total   = [t_ff_vec(i_ff_sel);   t_burn_shifted];
h_total   = [h_ff(i_ff_sel);       traj_burn.h_d];
vz_total  = [vz_ff(i_ff_sel);      traj_burn.vz_d];
az_total  = [az_ff(i_ff_sel);      traj_burn.az_d];

% Lateral X.
x_total   = [x_ff(i_ff_sel);       traj_burn.x_d];
vx_total  = [vx_ff(i_ff_sel);      traj_burn.vx_d];
ax_total  = [ax_ff(i_ff_sel);      traj_burn.ax_d];

% Lateral Y.
y_total   = [y_ff(i_ff_sel);       traj_burn.y_d];
vy_total  = [vy_ff(i_ff_sel);      traj_burn.vy_d];
ay_total  = [ay_ff(i_ff_sel);      traj_burn.ay_d];

% Atitude e empuxo.
theta_total = [theta_ff(i_ff_sel);  traj_burn.theta_d];
phi_total   = [phi_ff(i_ff_sel);    traj_burn.phi_d];
Tff_total   = [Tff_ff(i_ff_sel);    traj_burn.Tff_d];

% =========================================================================
%  CONVERSAO DAS VELOCIDADES DE REFERENCIA PARA O REFERENCIAL DO CORPO
% -------------------------------------------------------------------------
%  A malha de controle opera com velocidades no referencial do corpo
%  (u, v, w). Por isso convertemos as velocidades inerciais de referencia
%  de volta para o corpo. Como a curva de pouso so calcula phi e theta, assume-se psi nula durante o burn.
% =========================================================================
psi_total = [drop.psi(i_ff_sel); zeros(length(traj_burn.time), 1)];

% Prealoca os vetores de velocidade no referencial do corpo.
u_total = zeros(size(vx_total));
v_total = zeros(size(vy_total));
w_total = zeros(size(vz_total));

for i = 1:length(t_total)
    phi_i   = phi_total(i);
    th_i    = theta_total(i);
    psi_i   = psi_total(i);

    % Matriz R_IB (corpo -> inercial), identica a usada no drop test.
    R_IB = [cos(th_i)*cos(psi_i), sin(phi_i)*sin(th_i)*cos(psi_i) - cos(phi_i)*sin(psi_i), cos(phi_i)*sin(th_i)*cos(psi_i) + sin(phi_i)*sin(psi_i);
            cos(th_i)*sin(psi_i), sin(phi_i)*sin(th_i)*sin(psi_i) + cos(phi_i)*cos(psi_i), cos(phi_i)*sin(th_i)*sin(psi_i) - sin(phi_i)*cos(psi_i);
           -sin(th_i),            sin(phi_i)*cos(th_i),                                  cos(phi_i)*cos(th_i)];

    % R_BI (inercial -> corpo) e a transposta de R_IB (matriz ortogonal).
    R_BI = R_IB';

    V_I = [vx_total(i); vy_total(i); vz_total(i)];
    V_B = R_BI * V_I;

    u_total(i) = V_B(1);
    v_total(i) = V_B(2);
    w_total(i) = V_B(3);
end

% =========================================================================
%  MONTAGEM DA STRUCT DE SAIDA 'traj'
% -------------------------------------------------------------------------
%  Reune todos os perfis de referencia (vetores e timeseries) em uma unica
%  struct, pronta para ser consumida pelo modelo de controle no Simulink.
% =========================================================================

traj.T       = t_ff + T_burn;     % duracao total (queda + burn)
traj.T_burn  = T_burn;            % duracao apenas do burn
traj.t_ign   = t_ign;             % instante de ignicao
traj.h_ign   = h_ign;             % altitude de ignicao
traj.v_ign   = vz_ign;            % velocidade vertical na ignicao
traj.dt      = dt;
traj.time    = t_total;

traj.coef    = traj_burn.coef;    % coeficientes do polinomio vertical do burn
traj.s       = traj_burn.s;       % variavel normalizada do burn

% Perfis verticais (Z / altitude) e respectivas timeseries.
traj.h_d     = h_total;
traj.vz_d    = vz_total;
traj.az_d    = az_total;
traj.h_ref   = timeseries(h_total,  t_total, 'Name', 'h_ref');
traj.vz_ref  = timeseries(vz_total, t_total, 'Name', 'vz_ref');
traj.az_ref  = timeseries(az_total, t_total, 'Name', 'az_ref');

% Perfis laterais X.
traj.coef_x  = traj_burn.coef_x;
traj.x_d     = x_total;
traj.vx_d    = vx_total;
traj.ax_d    = ax_total;
traj.x_ref   = timeseries(x_total,  t_total, 'Name', 'x_ref');
traj.vx_ref  = timeseries(vx_total, t_total, 'Name', 'vx_ref');
traj.ax_ref  = timeseries(ax_total, t_total, 'Name', 'ax_ref');

% Perfis laterais Y.
traj.coef_y  = traj_burn.coef_y;
traj.y_d     = y_total;
traj.vy_d    = vy_total;
traj.ay_d    = ay_total;
traj.y_ref   = timeseries(y_total,  t_total, 'Name', 'y_ref');
traj.vy_ref  = timeseries(vy_total, t_total, 'Name', 'vy_ref');
traj.ay_ref  = timeseries(ay_total, t_total, 'Name', 'ay_ref');

% Perfis de atitude e empuxo de referencia.
traj.theta_d   = theta_total;
traj.phi_d     = phi_total;
traj.Tff_d     = Tff_total;
traj.theta_ref = timeseries(theta_total, t_total, 'Name', 'theta_ref');
traj.phi_ref   = timeseries(phi_total,   t_total, 'Name', 'phi_ref');
traj.Tff_ref   = timeseries(Tff_total,   t_total, 'Name', 'Tff_ref');

% Velocidades de referencia no referencial do corpo.
traj.u_d    = u_total;
traj.v_d    = v_total;
traj.w_d    = w_total;
traj.u_ref  = timeseries(u_total, t_total, 'Name', 'u_ref');
traj.v_ref  = timeseries(v_total, t_total, 'Name', 'v_ref');
traj.w_ref  = timeseries(w_total, t_total, 'Name', 'w_ref');

% ── Anexa os dados do envelope (drop test + scan) a struct de saida ──────
% Mantem disponivel para diagnostico e para os graficos posteriores.
envelope.drop          = drop;
envelope.scan_t        = scan_t;
envelope.scan_h        = scan_h;
envelope.scan_h_min    = scan_h_min;
envelope.scan_T_viol   = scan_T_viol;
envelope.scan_valid    = scan_valid;
envelope.t_early       = t_early;
envelope.t_late        = t_late;
envelope.h_early       = h_early;
envelope.h_late        = h_late;
envelope.alpha_ign     = alpha_ign;

traj.envelope = envelope;

% ── Diagnostico do empuxo exigido durante o burn ─────────────────────────
% Calcula o empuxo necessario ponto a ponto e a fracao do burn em que ele
% ficaria abaixo do empuxo minimo do motor (situacao indesejada).
T_req_burn = m0_pouso * sqrt(traj_burn.ax_d.^2 + traj_burn.ay_d.^2 ...
             + (g + traj_burn.az_d).^2);
frac_below = sum(T_req_burn < T_min_motor) / length(T_req_burn) * 100;
traj.frac_T_below_Tmin = frac_below;

if frac_below > 0
    fprintf('    [AVISO] %.1f%% do burn com T_req < T_min\n', frac_below);
else
    fprintf('    T_req >= T_min em 100%% do burn ✓\n');
end

% =========================================================================
%  DIAGNOSTICO FINAL (impressao do resumo da trajetoria gerada)
% =========================================================================

fprintf('\nTrajetória gerada (v6 — True-State):\n');
fprintf('  T_total = %.2f s  (queda=%.2f + burn=%.2f)\n', traj.T, t_ff, T_burn);
fprintf('  Pontos:           %d\n',        length(traj.time));
fprintf('  Pico az (burn):   %.2f m/s²\n', max(traj_burn.az_d));
fprintf('  Pico |theta_ref|: %.2f°\n',     rad2deg(max(abs(traj.theta_d))));
fprintf('  Pico |phi_ref|:   %.2f°\n',     rad2deg(max(abs(traj.phi_d))));
fprintf('  Pico T_ff (burn): %.1f kN  (T_max=%.0f kN)\n', ...
        max(traj_burn.Tff_d)/1e3, empuxo_max/1e3);
fprintf('  Posição final:    (%.4f, %.4f, %.4f) m\n', ...
        traj.x_d(end), traj.y_d(end), traj.h_d(end));
fprintf('  Velocidade final: (%.4f, %.4f, %.4f) m/s\n', ...
        traj.vx_d(end), traj.vy_d(end), traj.vz_d(end));

% ── Validacao das condicoes de contorno do burn ──────────────────────────
% Verifica se a trajetoria realmente termina nas condicoes de pouso
% desejadas (as CI globais nao se aplicam aqui, pois vem da simulacao).
fprintf('\nValidação das condições de contorno (burn):\n');
tol = 1e-3;   % tolerancia relaxada (CI vem da simulacao, nao de formula)
check_bc('h(T)',   traj.h_d(end),  hf,  tol);
check_bc('vz(T)',  traj.vz_d(end), vzf, tol);
check_bc('x(T)',   traj.x_d(end),  xf,  tol);
check_bc('y(T)',   traj.y_d(end),  yf,  tol);

% =========================================================================
%  GRAFICOS — Perfis temporais da trajetoria
% -------------------------------------------------------------------------
%  Painel 5x3 com posicao, velocidades (corpo e inercial), aceleracoes,
%  atitude e empuxo ao longo do tempo. A linha verde tracejada marca o
%  instante de ignicao em todos os subgraficos.
% =========================================================================
az_max_fis = (empuxo_max / m0_pouso) - g;   % aceleracao vertical maxima fisica
az_alvo    = margem * az_max_fis;           % aceleracao vertical alvo (com margem)

% Janela ampla para acomodar confortavelmente as 5 linhas de graficos.
figure('Name','Guidance - 3D Trajectory', ...
       'NumberTitle','off','Position',[40 40 1400 900]);

% ── LINHA 1: POSICAO INERCIAL (x, y, h) ──────────────────────────────────
subplot(5,3,1);
plot(traj.time, traj.x_d, 'b-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff, 'g--', 'ignition', 'LineWidth', 1);
ylabel('x_d [m]'); grid on; title('X (Lateral)');

subplot(5,3,2);
plot(traj.time, traj.y_d, 'b-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff, 'g--', 'LineWidth', 1);
ylabel('y_d [m]'); grid on; title('Y (Lateral)');

subplot(5,3,3);
plot(traj.time, traj.h_d, 'b-', 'LineWidth', 1.8); hold on;
yline(0,'k--','Ground','LineWidth',0.8); xline(t_ff, 'g--', 'ignition', 'LineWidth', 1);
ylabel('h_d [m]'); grid on; title('Z / h (Vertical)');

% ── LINHA 2: VELOCIDADE NO CORPO (u, v, w) ───────────────────────────────
subplot(5,3,4);
plot(traj.time, traj.u_d, 'Color', [0 0.6 0.8], 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
ylabel('u_d (Body) [m/s]'); grid on; title('Body Velocity');

subplot(5,3,5);
plot(traj.time, traj.v_d, 'Color', [0 0.6 0.8], 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
ylabel('v_d (Body) [m/s]'); grid on; title('Body Velocity');

subplot(5,3,6);
plot(traj.time, traj.w_d, 'Color', [0 0.6 0.8], 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','ignition','LineWidth',1);
ylabel('w_d (Body) [m/s]'); grid on; title('Body Velocity');

% ── LINHA 3: VELOCIDADE INERCIAL (vx, vy, vz) ────────────────────────────
subplot(5,3,7);
plot(traj.time, traj.vx_d, 'r-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
ylabel('vx_d (Inerc) [m/s]'); grid on; title('Inercial Velocity');

subplot(5,3,8);
plot(traj.time, traj.vy_d, 'r-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
ylabel('vy_d (Inerc) [m/s]'); grid on; title('Inercial Velocity');

subplot(5,3,9);
plot(traj.time, traj.vz_d, 'r-', 'LineWidth', 1.8); hold on;
yline(0,'k--','v=0','LineWidth',0.8); xline(t_ff,'g--','ignition','LineWidth',1);
ylabel('vz_d (Inerc) [m/s]'); grid on; title('Inercial Velocity');

% ── LINHA 4: ACELERACAO INERCIAL (ax, ay, az) ────────────────────────────
subplot(5,3,10);
plot(traj.time, traj.ax_d, 'g-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
ylabel('ax_d [m/s^2]'); grid on; title('Aceleration');

subplot(5,3,11);
plot(traj.time, traj.ay_d, 'g-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
ylabel('ay_d [m/s^2]'); grid on; title('Aceleration');

subplot(5,3,12);
plot(traj.time, traj.az_d, 'g-', 'LineWidth', 1.8); hold on;
% Linhas de referencia: aceleracao maxima fisica e aceleracao alvo.
yline(az_max_fis, 'r--', sprintf('az_{max}=%.1f',az_max_fis), 'LineWidth', 1.2);
yl = yline(az_alvo, '--', sprintf('az_{objctive}=%.1f',az_alvo), 'LineWidth', 1.2);
yl.Color = [0.7 0.7 0];
yline(0,'k--','LineWidth',0.6); xline(t_ff,'g--','ignition','LineWidth',1);
ylabel('az_d [m/s^2]'); grid on; title('Aceleration');

% ── LINHA 5: ATITUDE E EMPUXO ────────────────────────────────────────────
subplot(5,3,13);
plot(traj.time, rad2deg(traj.theta_d), 'm-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
xlabel('Time [s]'); ylabel('\theta_{ref} [deg]'); grid on;
title('Yaw (ref atitude)');

subplot(5,3,14);
plot(traj.time, rad2deg(traj.phi_d), 'm-', 'LineWidth', 1.8); hold on;
yline(0,'k--','LineWidth',0.8); xline(t_ff,'g--','LineWidth',1);
xlabel('Time [s]'); ylabel('\phi_{ref} [deg]'); grid on;
title('Pitch (ref atitude)');

subplot(5,3,15);
plot(traj.time, traj.Tff_d/1e3, 'k-', 'LineWidth', 1.8); hold on;
% Limites superior (empuxo maximo) e inferior (empuxo minimo) do motor.
yline(empuxo_max/1e3, 'r--', sprintf('T_{max}=%.0fkN',empuxo_max/1e3), 'LineWidth', 1.2);
yline(T_min_motor/1e3, 'b--', sprintf('T_{min}=%.0fkN',T_min_motor/1e3),'LineWidth',1.2);
xline(t_ff,'g--','ignition','LineWidth',1);
xlabel('Time [s]'); ylabel('T_{ff} [kN]'); grid on;
title('Feedforward Thrust');

sgtitle(sprintf(['Suicide Burn 3D  |  t_{ign}=%.2fs  h_{ign}=%.0fm  |  ' ...
         'T_{burn}=%.1fs  |  \\alpha=%.2f  [%.2f–%.2f s]'], ...
         t_ign, h_ign, T_burn, alpha_ign, t_early, t_late), ...
        'FontSize',12,'FontWeight','bold');

% =========================================================================
%  GRAFICO — Trajetoria 3D no espaco
% -------------------------------------------------------------------------
%  Mostra o caminho completo do foguete no espaco (x, y, h), o ponto de
%  ignicao, o inicio da queda, o alvo de pouso e a sombra projetada no solo.
% =========================================================================
figure('Name','3D Trajectory','NumberTitle','off',...
       'Position',[100 100 800 650]);

% Linha grossa para a trajetoria; 'axis equal' foi propositalmente omitido
% para nao distorcer a escala vertical e prejudicar o render da camera.
plot3(traj.x_d, traj.y_d, traj.h_d, 'b-', 'LineWidth', 3.5); hold on;

% Marcadores destacando os pontos notaveis da trajetoria.
scatter3(x_ign, y_ign, h_ign, 120, 'g', 'filled', 'd', 'MarkerEdgeColor', 'k');
scatter3(x0, y0, h0, 90, 'b', 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter3(xf, yf, hf, 180, 'r', 'filled', 'p', 'MarkerEdgeColor', 'k');

% Sombra da trajetoria projetada sobre o solo (h = 0).
plot3(traj.x_d, traj.y_d, zeros(size(traj.h_d)), 'Color', [0.4 0.4 0.4], 'LineStyle', '--', 'LineWidth', 1.5);

grid on;
% axis equal;  % intencionalmente removido para nao "esmagar" o grafico
view(45, 25);
xlabel('X Axis [m]', 'FontWeight', 'bold');
ylabel('Y Axis [m]', 'FontWeight', 'bold');
zlabel('Altitude h [m]', 'FontWeight', 'bold');
legend('3D Trajectory','Ignition Point','Free-Fall start','Landing objective (0,0,0)','Ground shadow', ...
       'Location','northeast', 'FontSize', 10);
title(sprintf('3D Spatial Reference  |  Ignition Altitude: %.0f m', h_ign), 'FontSize', 12);

% =========================================================================
%  GRAFICO — Drop test e True-State Scan
% -------------------------------------------------------------------------
%  Quatro paineis de diagnostico: altitude da queda com a janela de
%  ignicao, desvio lateral real vs analitico, e os dois criterios do scan
%  (altitude minima do polinomio e violacao de empuxo minimo).
% =========================================================================
figure('Name','Drop Test & True-State Scan','NumberTitle','off',...
       'Position',[180 180 1100 500]);

% ── Painel 1: altitude da queda e janela de ignicao ─────────────────────
subplot(2,2,1);
plot(drop.t, drop.h, 'b-', 'LineWidth', 1.8); hold on;
yline(0, 'k--', 'Ground', 'LineWidth', 0.8);
xline(t_early, 'g--', sprintf('t_{early}=%.2fs', t_early), 'LineWidth', 1.2);
xline(t_late,  'r--', sprintf('t_{late}=%.2fs', t_late),   'LineWidth', 1.2);
xline(t_ign,   'k-',  sprintf('t_{ign}=%.2fs', t_ign),     'LineWidth', 1.8);
% Sombreia a regiao correspondente a janela de ignicao viavel.
yl = ylim;
patch([t_early t_late t_late t_early], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.2 0.8 0.2], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
xlabel('Time [s]'); ylabel('h [m]'); grid on;
title('Drop Test: Altitude');
legend('h(t)', 'Ground', 't_{early}', 't_{late}', 't_{ign}', 'Window', ...
       'Location','northeast');

% ── Painel 2: desvio lateral real vs analitico ──────────────────────────
subplot(2,2,2);
plot(drop.t, drop.x - x0, 'b-', 'LineWidth', 1.5); hold on;
plot(drop.t, drop.y - y0, 'r-', 'LineWidth', 1.5);
% Referencia analitica (movimento retilineo uniforme) para comparacao.
plot(drop.t, vx0 * drop.t, 'b--', 'LineWidth', 0.8);
plot(drop.t, vy0 * drop.t, 'r--', 'LineWidth', 0.8);
xline(t_ign, 'k-', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Δposition [m]'); grid on;
title('Lateral Drift: real vs analitic');
legend('Δx real', 'Δy real', 'Δx analitic', 'Δy analitic', 't_{ign}', ...
       'Location','northwest');

% ── Painel 3: altitude minima do polinomio (verifica furo no solo) ──────
subplot(2,2,3);
plot(scan_t, scan_h_min, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3); hold on;
yline(h_floor, 'r--', sprintf('h_{floor}=%d m', h_floor), 'LineWidth', 1.2);
xline(t_early, 'g--', 'LineWidth', 1.2);
xline(t_late,  'r--', 'LineWidth', 1.2);
xline(t_ign,   'k-',  'LineWidth', 1.8);
xlabel('t_{ign} candidate [s]'); ylabel('polynomial min(h_d) [m]'); grid on;
title('Trajectory-undershoot check');

% ── Painel 4: violacao de empuxo minimo ao longo do scan ────────────────
subplot(2,2,4);
plot(scan_t, scan_T_viol*100, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 3); hold on;
yline(T_VIOL_TOL*100, 'b--', sprintf('Tol=%.0f%%', T_VIOL_TOL*100), 'LineWidth', 1.2);
xline(t_early, 'g--', 'LineWidth', 1.2);
xline(t_late,  'r--', 'LineWidth', 1.2);
xline(t_ign,   'k-',  'LineWidth', 1.8);
xlabel('t_{ign} candidate [s]'); ylabel('T_{viol} [% do burn]'); grid on;
title('T_{min} Violation');

sgtitle('Drop Test & True-State Scan', 'FontSize', 11, 'FontWeight', 'bold');

fprintf('\nPronto. ''traj'', ''envelope'', ''drop'' disponíveis no workspace.\n');
fprintf('Campos: traj.envelope.drop, .t_early, .t_late, .scan_*\n\n');

% =========================================================================
% =========================================================================
%  FUNCOES LOCAIS
%  Definidas ao final do script (escopo restrito a este arquivo).
% =========================================================================
% =========================================================================

% -------------------------------------------------------------------------
%  scan_frame_rapido
%  Avalia um unico frame candidato de ignicao (usado dentro do scan).
%
%  Recebe o estado REAL do foguete como condicao inicial e tenta gerar a
%  curva de pouso 3D acoplada. Para ganhar desempenho, usa avaliacao em
%  curto-circuito (short-circuit):
%    1. Resolve primeiro o canal vertical (Z).
%    2. Se a altitude minima ja for negativa (fura o solo), retorna sem
%       resolver os canais laterais.
%    3. So calcula X, Y, empuxo 3D e violacao de empuxo minimo se Z passar.
%
%  ENTRADAS (CI e condicoes finais nos tres eixos, mais parametros fisicos).
%  SAIDAS:
%    ok     - true se o calculo do tempo de queima (fzero) convergiu.
%    h_min  - altitude minima do polinomio de burn.
%    T_viol - fracao [0,1] do burn com empuxo exigido abaixo do minimo.
% -------------------------------------------------------------------------
function [ok, h_min, T_viol] = scan_frame_rapido( ...
    h0,vz0,az0, hf,vzf,azf, ...
    x0,vx0,ax0, xf,vxf,axf, ...
    y0,vy0,ay0, yf,vyf,ayf, ...
    empuxo_max, T_min_motor, m0, margem, dt, g)

    % Valores de retorno padrao (caso algo falhe).
    ok     = false;
    h_min  = NaN;
    T_viol = NaN;

    % Distancia vertical a percorrer; se nao for positiva, nao ha manobra.
    dh = h0 - hf;
    if dh <= 0
        return;
    end

    % ── Encontra o tempo de queima T_ot via empuxo 3D acoplado ──────────
    % Busca o T cujo pico de empuxo iguala a margem do empuxo maximo.
    % Tolerancia relaxada (1e-3) por se tratar de uma triagem rapida.
    try
        f_res = @(T_) max_thrust_dado_T(T_, ...
            h0,vz0,az0, hf,vzf,azf, dh, dt, ...
            x0,vx0,ax0, xf,vxf,axf, ...
            y0,vy0,ay0, yf,vyf,ayf, ...
            m0, g) - (margem * empuxo_max);

        opts = optimset('TolX',1e-3,'Display','off');
        T_ot = fzero(f_res, [5, 300], opts);
    catch
        return;   % nao convergiu: frame invalido
    end

    % ── Curto-circuito: resolve primeiro o canal vertical Z ─────────────
    try
        [~, t, ~, ~, ~, ~, h_d, ~, az_d] = ...
            resolver_scurve_z(T_ot, h0,vz0,az0,hf,vzf,azf,dh,dt);
    catch
        return;
    end

    h_min = min(h_d);
    ok    = true;

    % Se a trajetoria fura o solo, marca como violacao total e retorna sem
    % gastar tempo resolvendo os canais laterais.
    if h_min < 0
        T_viol = 1.0;     % 100% de violacao (frame reprovado)
        return;
    end

    % ── Canais laterais e empuxo 3D (so se o canal Z passou) ────────────
    try
        [~, ~, ~, ax_d] = resolver_canal(T_ot, x0,vx0,ax0, xf,vxf,axf, t);
        [~, ~, ~, ay_d] = resolver_canal(T_ot, y0,vy0,ay0, yf,vyf,ayf, t);
    catch
        T_viol = 1.0;
        return;
    end

    % Modulo do empuxo necessario ponto a ponto (norma do vetor de forca).
    Tff_d = m0 * sqrt(ax_d.^2 + ay_d.^2 + (g + az_d).^2);

    % ── Fracao do burn em que o empuxo exigido fica abaixo do minimo ────
    T_viol = sum(Tff_d < T_min_motor) / length(Tff_d);
end

% -------------------------------------------------------------------------
%  scurve_calc_3d
%  Calcula a trajetoria de pouso completa do burn, com os tres eixos
%  acoplados pelo empuxo total.
%
%  Otimiza o tempo de queima T de modo que o pico de empuxo necessario
%  iguale a margem do empuxo maximo, e entao gera os polinomios de
%  referencia para Z, X e Y, alem da atitude (theta, phi) e do empuxo de
%  feedforward (Tff). Retorna tudo na struct 'traj'.
% -------------------------------------------------------------------------
function traj = scurve_calc_3d( ...
    h0,vz0,az0, hf,vzf,azf, ...
    x0,vx0,ax0, xf,vxf,axf, ...
    y0,vy0,ay0, yf,vyf,ayf, ...
    empuxo_max, T_min_motor, m0, margem, dt, g) %#ok<INUSD>

    dh = h0 - hf;   % distancia vertical a percorrer no burn

    % Funcao-residuo: diferenca entre o pico de empuxo (para um dado T) e
    % a margem do empuxo maximo. A raiz fornece o T otimo.
    f_res = @(T_) max_thrust_dado_T(T_, ...
        h0,vz0,az0, hf,vzf,azf, dh, dt, ...
        x0,vx0,ax0, xf,vxf,axf, ...
        y0,vy0,ay0, yf,vyf,ayf, ...
        m0, g) - (margem * empuxo_max);

    opts = optimset('TolX',1e-4,'Display','off');   % tolerancia mais fina aqui
    T_ot = fzero(f_res, [5, 300], opts);            % busca T no intervalo [5, 300] s

    % Gera o perfil vertical (Z) com o T otimo.
    [coef_z_norm, t, ~, s, ~, ~, h_d, vz_d, az_d] = ...
        resolver_scurve_z(T_ot, h0,vz0,az0,hf,vzf,azf,dh,dt);

    % Gera os perfis laterais X e Y sobre a mesma base de tempo.
    [coef_x, x_d, vx_d, ax_d] = resolver_canal( ...
        T_ot, x0, vx0, ax0, xf, vxf, axf, t);

    [coef_y, y_d, vy_d, ay_d] = resolver_canal( ...
        T_ot, y0, vy0, ay0, yf, vyf, ayf, t);

    % Atitude de referencia obtida da direcao do vetor de aceleracao.
    theta_d = atan2(ax_d, g + az_d);
    phi_d   = -atan2(ay_d, g + az_d);
    % Empuxo de feedforward (norma do vetor de forca exigido).
    Tff_d   = m0 * sqrt(ax_d.^2 + ay_d.^2 + (g + az_d).^2);

    % Empacota vetores numericos.
    traj.T       = T_ot;
    traj.dt      = dt;
    traj.time    = t;
    traj.coef    = coef_z_norm;
    traj.s       = s;
    traj.h_d     = h_d;      traj.vz_d  = vz_d;    traj.az_d  = az_d;
    traj.coef_x  = coef_x;
    traj.x_d     = x_d;      traj.vx_d  = vx_d;    traj.ax_d  = ax_d;
    traj.coef_y  = coef_y;
    traj.y_d     = y_d;      traj.vy_d  = vy_d;    traj.ay_d  = ay_d;
    traj.theta_d = theta_d;  traj.phi_d = phi_d;   traj.Tff_d = Tff_d;

    % Empacota as mesmas grandezas como timeseries (uso direto no Simulink).
    traj.h_ref     = timeseries(h_d,     t, 'Name','h_ref');
    traj.vz_ref    = timeseries(vz_d,    t, 'Name','vz_ref');
    traj.az_ref    = timeseries(az_d,    t, 'Name','az_ref');
    traj.x_ref     = timeseries(x_d,     t, 'Name','x_ref');
    traj.vx_ref    = timeseries(vx_d,    t, 'Name','vx_ref');
    traj.ax_ref    = timeseries(ax_d,    t, 'Name','ax_ref');
    traj.y_ref     = timeseries(y_d,     t, 'Name','y_ref');
    traj.vy_ref    = timeseries(vy_d,    t, 'Name','vy_ref');
    traj.ay_ref    = timeseries(ay_d,    t, 'Name','ay_ref');
    traj.theta_ref = timeseries(theta_d, t, 'Name','theta_ref');
    traj.phi_ref   = timeseries(phi_d,   t, 'Name','phi_ref');
    traj.Tff_ref   = timeseries(Tff_d,   t, 'Name','Tff_ref');
end

% -------------------------------------------------------------------------
%  max_thrust_dado_T
%  Para um tempo de queima T candidato, calcula o pico do empuxo 3D
%  necessario ao longo de todo o burn. E a funcao-objetivo usada pelo
%  fzero na otimizacao de T (tanto no scan quanto na geracao final).
%  Retorna NaN se algum dos canais nao puder ser resolvido.
% -------------------------------------------------------------------------
function T_pk = max_thrust_dado_T(T_, ...
    h0,vz0,az0, hf,vzf,azf, dh, dt, ...
    x0,vx0,ax0, xf,vxf,axf, ...
    y0,vy0,ay0, yf,vyf,ayf, ...
    m0, g)

    try
        [~, t, ~, ~, ~, ~, ~, ~, az_d] = ...
            resolver_scurve_z(T_, h0,vz0,az0,hf,vzf,azf,dh,dt);
        [~, ~, ~, ax_d] = resolver_canal(T_, x0,vx0,ax0, xf,vxf,axf, t);
        [~, ~, ~, ay_d] = resolver_canal(T_, y0,vy0,ay0, yf,vyf,ayf, t);
        Tff_d = m0 * sqrt(ax_d.^2 + ay_d.^2 + (g + az_d).^2);
        T_pk = max(Tff_d);
    catch
        T_pk = NaN;
    end
end

% -------------------------------------------------------------------------
%  az_pico_dado_T
%  Versao simplificada que retorna apenas o pico da aceleracao vertical
%  para um T candidato. Mantida para compatibilidade e depuracao; nao e
%  chamada no fluxo principal atual.
% -------------------------------------------------------------------------
function az_pk = az_pico_dado_T(T_, h0,vz0,az0,hf,vzf,azf,dh,dt) %#ok<DEFNU>
    try
        [~,~,~,~,~,~,~,~,az_d_] = ...
            resolver_scurve_z(T_,h0,vz0,az0,hf,vzf,azf,dh,dt);
        az_pk = max(az_d_);
    catch
        az_pk = NaN;
    end
end

% -------------------------------------------------------------------------
%  resolver_scurve_z
%  Monta e resolve o polinomio de 5a ordem do eixo vertical (Z).
%
%  Trabalha com uma variavel normalizada s in [0,1] e seu parametro tau =
%  t/T. As seis condicoes de contorno (posicao, velocidade e aceleracao no
%  inicio e no fim) determinam os seis coeficientes do polinomio via a
%  matriz M. A partir de s e suas derivadas, reconstroi h, vz e az reais.
% -------------------------------------------------------------------------
function [coef,t,tau,s,s_d1,s_d2,h_d,vz_d,az_d] = ...
         resolver_scurve_z(T_, h0,vz0,az0,hf,vzf,azf,dh,dt)

    % Condicoes de contorno na variavel normalizada s (inicio = 0, fim = 1).
    s0   =  0;
    s1   =  1;
    sd0  = -vz0 * T_  / dh;       % velocidade normalizada inicial
    sd1  = -vzf * T_  / dh;       % velocidade normalizada final
    sdd0 = -az0 * T_^2 / dh;      % aceleracao normalizada inicial
    sdd1 = -azf * T_^2 / dh;      % aceleracao normalizada final

    % Matriz que relaciona os coeficientes do polinomio as condicoes de
    % contorno (avaliadas em tau=0 e tau=1).
    M = [1  0  0   0   0    0  ;
         0  1  0   0   0    0  ;
         0  0  2   0   0    0  ;
         1  1  1   1   1    1  ;
         0  1  2   3   4    5  ;
         0  0  2   6   12   20 ];

    b    = [s0; sd0; sdd0; s1; sd1; sdd1];
    coef = M \ b;        % resolve o sistema linear para os coeficientes

    t   = (0:dt:T_)';    % base de tempo discreta do burn
    tau = t / T_;        % tempo normalizado em [0,1]

    % Polinomio normalizado s(tau) e suas duas primeiras derivadas.
    s    = coef(1) + coef(2)*tau      + coef(3)*tau.^2   + ...
           coef(4)*tau.^3 + coef(5)*tau.^4 + coef(6)*tau.^5;
    s_d1 =           coef(2)          + 2*coef(3)*tau    + ...
           3*coef(4)*tau.^2 + 4*coef(5)*tau.^3 + 5*coef(6)*tau.^4;
    s_d2 =                     2*coef(3)                 + ...
           6*coef(4)*tau    + 12*coef(5)*tau.^2 + 20*coef(6)*tau.^3;

    % Converte de volta para as grandezas fisicas reais (altitude,
    % velocidade e aceleracao verticais).
    h_d  =  h0 - dh * s;
    vz_d = -(dh / T_)   * s_d1;
    az_d = -(dh / T_^2) * s_d2;
end

% -------------------------------------------------------------------------
%  resolver_canal
%  Monta e resolve o polinomio de 5a ordem para um canal lateral (X ou Y).
%
%  Diferente do canal vertical, aqui o polinomio e expresso diretamente na
%  variavel fisica (posicao), com as condicoes de contorno escaladas por T
%  e T^2 para velocidade e aceleracao. Retorna posicao, velocidade e
%  aceleracao desejadas sobre a base de tempo fornecida.
% -------------------------------------------------------------------------
function [coef, p_d, v_d, a_d] = resolver_canal(T_, p0, v0, a0, pf, vf, af, t)

    tau = t / T_;        % tempo normalizado em [0,1]

    % Mesma matriz de condicoes de contorno do canal vertical.
    M = [1  0  0   0   0    0  ;
         0  1  0   0   0    0  ;
         0  0  2   0   0    0  ;
         1  1  1   1   1    1  ;
         0  1  2   3   4    5  ;
         0  0  2   6   12   20 ];

    % Vetor de condicoes de contorno: posicao, velocidade*T e aceleracao*T^2
    % no inicio e no fim.
    b    = [p0; T_*v0; T_^2*a0; pf; T_*vf; T_^2*af];
    coef = M \ b;

    % Polinomio de posicao e suas derivadas (reescaladas por T e T^2).
    p_d = coef(1) + coef(2)*tau      + coef(3)*tau.^2   + ...
          coef(4)*tau.^3 + coef(5)*tau.^4 + coef(6)*tau.^5;
    v_d = (coef(2) + 2*coef(3)*tau   + 3*coef(4)*tau.^2 + ...
           4*coef(5)*tau.^3 + 5*coef(6)*tau.^4) / T_;
    a_d = (2*coef(3) + 6*coef(4)*tau + 12*coef(5)*tau.^2 + ...
           20*coef(6)*tau.^3) / T_^2;
end

% -------------------------------------------------------------------------
%  Utilitarios
% -------------------------------------------------------------------------

% local_get
%   Le um campo de uma struct de forma segura: retorna o valor do campo se
%   ele existir, ou um valor padrao caso contrario. Evita erros quando a
%   struct de parametros nao define todos os campos opcionais.
function v = local_get(s, field, default)
    if isstruct(s) && isfield(s, field)
        v = s.(field);
    else
        v = default;
    end
end

% check_bc
%   Verifica se um valor calculado bate com o valor esperado dentro de uma
%   tolerancia e imprime "OK" ou "ERRO" de forma padronizada. Usada na
%   validacao das condicoes de contorno do burn.
function check_bc(nome, valor, esperado, tol)
    if abs(valor - esperado) <= tol
        fprintf('  OK  %-8s = %+9.4f  (esperado: %+.4f)\n', nome, valor, esperado);
    else
        fprintf('  ERRO %-7s = %+9.4f  (esperado: %+.4f)  diff: %.2e\n', ...
                nome, valor, esperado, abs(valor-esperado));
    end
end

% ev_solo
%   Funcao de evento para integradores (ode*): detecta o instante em que a
%   altitude (estado 3) cruza zero, descendo (impacto no solo), e encerra a
%   integracao. Mantida no arquivo para uso por simulacoes que a requeiram.
function [value, isterminal, direction] = ev_solo(~, x) %#ok<DEFNU>
    value      = x(3);     % altitude; o evento dispara quando cruza zero
    isterminal = 1;        % interrompe a integracao no evento
    direction  = -1;       % apenas na descida (valor decrescente)
end