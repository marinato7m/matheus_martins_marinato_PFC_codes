% =========================================================================
% OTIMIZAR_GANHOS_2ESTAGIOS.M
% -------------------------------------------------------------------------
%
% DESCRIÇÃO GERAL
%   Rotina de tuning automático dos pesos das matrizes Q e R do
%   controlador MPC de pouso do veículo. O problema de tuning é
%   tratado como uma otimização sobre 7 multiplicadores adimensionais
%   ("knobs"), cada um escalando um grupo de termos das matrizes de
%   ponderação a partir de uma referência inicial de Bryson.
% MOTIVAÇÃO DO PIPELINE EM 2 ESTÁGIOS
%   Bayesian Optimization (BO) é eficiente em poucas dimensões, mas degrada
%   em espaços 7D com bounds largos: o modelo gaussiano gasta avaliações
%   apenas mapeando regiões irrelevantes. A estratégia aqui é primeiro
%   varrer o espaço completo de forma barata (LHS), descobrir quais knobs
%   realmente importam e onde estão os bons valores, e só então rodar o BO
%   num subespaço pequeno e bem delimitado, onde ele converge rápido.
%
% ─────────────────────────────────────────────────────────────────────────
% ESTÁGIO 1 — EXPLORAÇÃO (Latin Hypercube Sampling)
%   Gera 200 pontos bem distribuídos no espaço 7D via LHS.
%   Avalia todos, classifica, e faz análise de sensibilidade:
%     - Correlação de Spearman (quais knobs importam?)
%     - Análise de cluster dos TOP-30 (quais knobs já convergiram?)
%     - Gráficos de dependência parcial
%   Resultado: knobs estáveis são FIXADOS, bounds dos restantes APERTADOS.
%
%   A ideia central: se entre as melhores soluções um knob assume sempre
%   praticamente o mesmo valor (baixa dispersão), ele já "convergiu" e pode
%   ser congelado; se varia muito, ele ainda precisa ser otimizado.
%
% ESTÁGIO 2 — REFINAMENTO (Bayesian Optimization)
%   Opera em espaço reduzido com bounds estreitos.
%   O melhor ponto do LHS é injetado como chute inicial do BO, garantindo
%   que o refinamento parta de uma solução já razoável.
%
% ─────────────────────────────────────────────────────────────────────────
% FUNÇÃO-OBJETIVO (CUSTO)
%   Cada candidato é avaliado por um "mini Monte-Carlo": a mesma sintonia é
%   simulada sobre um conjunto fixo de condições iniciais de desafio (CIs),
%   e o custo agrega, via percentil 95, as métricas de pouso (erro lateral,
%   velocidade vertical, atitude e ângulo máximo), penalizando ainda os
%   casos que não pousam. Menor custo = melhor sintonia.
%
% PRÉ-REQUISITOS:
%   >> init_mpc            % carrega p, opts_qp
%   (gerar_traj_rapido.m e rotation_matrix.m no path)
%
% EXECUÇÃO:
%   >> init_mpc
%   >> otimizar_ganhos_2estagios
%
% SAÍDAS GERADAS:
%   resultados_estagio1.mat  — dados do LHS + análise de sensibilidade
%   resultados_2estagios.mat — resultado final consolidado (Q_opt, R_opt)
%   Bloco de texto pronto para colar no init_mpc.m
%
% TEMPO ESTIMADO: ~35-75 min total
%
% Matheus Martins Marinato — TCC UFSM 2025/2026
% =========================================================================

%% ═══════════════════════════════════════════════════════════════════════
%  0. VERIFICAÇÕES E SETUP
% ═══════════════════════════════════════════════════════════════════════
%   Confere se o ambiente foi corretamente inicializado por init_mpc antes
%   de prosseguir. Falhar cedo aqui evita erros lá na frente, já
%   que p e opts_qp são consumidos por praticamente todas as etapas.

% 'p' agrupa os parâmetros físicos/limites do veículo (massa, empuxo, TVC).
if ~exist('p', 'var')
    error('Variável "p" não encontrada. Execute init_mpc primeiro.');
end
% 'opts_qp' contém as opções do solver QP usado pelo MPC dentro do modelo.
if ~exist('opts_qp', 'var')
    error('Variável "opts_qp" não encontrada. Execute init_mpc primeiro.');
end
% A geração das trajetórias de referência depende desta função auxiliar.
if exist('gerar_traj_rapido', 'file') ~= 2
    error('gerar_traj_rapido.m não encontrado no path.');
end

fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║    OTIMIZAÇÃO 2-ESTÁGIOS: LHS + BAYESOPT (7 knobs, mini-MC) ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% ═══════════════════════════════════════════════════════════════════════
%  1. CONFIGURAÇÃO GLOBAL
% ═══════════════════════════════════════════════════════════════════════
%   Toda a parametrização do estudo fica concentrada na struct 'cfg', que é
%   passada adiante para as funções de avaliação. Centralizar aqui mantém o
%   experimento reprodutível e fácil de ajustar sem caçar constantes soltas
%   pelo corpo do script.

% Modelo Simulink a ser simulado. As duas primeiras linhas ficam comentadas
% como histórico das versões anteriores do modelo (sem navegação / etc.);
% a versão ativa inclui a dinâmica dos atuadores.
cfg.modelo      = 'rocket_6dof_sim';

% Quando true, as N_CI simulações de cada candidato rodam em paralelo
% (parsim), o que reduz bastante o tempo total da otimização.
cfg.USAR_PARSIM = true;

% ── Budget ──────────────────────────────────────────────────────────────
%   Orçamento de avaliações de cada estágio. É o principal fator do tempo
%   de execução, já que cada avaliação dispara N_CI simulações completas.
cfg.N_LHS       = 200;     % pontos do Latin Hypercube (estágio 1)
cfg.N_TOP       = 30;      % melhores pontos pra análise de cluster
cfg.N_BAYESOPT  = 50;      % avaliações do bayesopt (estágio 2)
cfg.N_SEED_BO   = 10;      % seed points do bayesopt

% ── Requisitos ──────────────────────────────────────────────────────────
%   Requisitos de desempenho de pouso (p95). Servem de referência para o
%   relatório final (marcação OK/FALHA); não entram diretamente no custo.
cfg.REQ_LAT       = 5.0;   % erro lateral máximo no toque [m]
cfg.REQ_VZ        = 3.0;   % velocidade vertical máxima no toque [m/s]
cfg.REQ_ATT       = 6.0;   % erro de atitude no toque [graus]
cfg.REQ_THETA_MAX = 20.0;  % pico de theta ao longo da descida [graus]

% ── Referências Bryson ──────────────────────────────────────────────────
%   Valores de "máximo aceitável" de cada estado, usados na regra de Bryson
%   para normalizar Q (peso ~ 1/valor_máximo^2). Os knobs apenas escalam
%   essas referências, de modo que a sintonia parte de uma base fisicamente
%   significativa em vez de números arbitrários.
cfg.bryson.pos  = 5.0;            % posição lateral [m]
cfg.bryson.vel  = 3.0;            % velocidade horizontal [m/s]
cfg.bryson.h    = 2.0;            % altitude (canal vertical de posição) [m]
cfg.bryson.vz   = 1.0;            % velocidade vertical [m/s]
cfg.bryson.ang  = deg2rad(5);     % ângulos de atitude [rad]
cfg.bryson.rate = deg2rad(5);     % taxas angulares [rad/s]
% Bases de R derivadas da faixa física dos atuadores (empuxo e TVC):
cfg.R_base_T    = 1 / (p.T_max - p.T_min)^2;   % normalização do empuxo
cfg.R_base_tvc  = 1 / (p.delta_max)^2;         % normalização do TVC
cfg.p           = p;   % guarda os parâmetros do veículo dentro de cfg

% ── Bounds dos 7 knobs (log-space) ─────────────────────────────────────
%   Limites de busca de cada multiplicador. Os knobs variam por ordens de
%   grandeza, então a amostragem e a otimização são feitas em escala log,
%   garantindo cobertura uniforme.
%   Cada linha: [nome, lb, ub]
cfg.knob_names = {'k_pos','k_vel','k_h','k_vz','k_ang','k_rate','k_tvc'};
cfg.N_KNOBS    = 7;
cfg.lb = [0.5,  0.1,  1,   1,   0.1,  0.05,  0.1];   % lower bounds
cfg.ub = [100,  50,   50,  100, 50,   25,    10 ];    % upper bounds

% ── Limiar de convergência ──────────────────────────────────────────────
%   Se o coef. de variação (std/mean) do knob nos TOP-N for menor que
%   este valor, o knob é considerado "convergido" e fixado no estágio 2.
%   Valor alto = conservador, fixa menos knobs (mais segurança, BO maior).
cfg.CV_THRESHOLD = 0.40;   % 40% — conservador (não fixa cedo demais)

%% ═══════════════════════════════════════════════════════════════════════
%  2. CONDIÇÕES INICIAIS DE DESAFIO
% ═══════════════════════════════════════════════════════════════════════
%   "mini Monte-Carlo". Premia robustez e evita sintonias que só funcionam num cenário.
%   Colunas: altitude, posições e velocidades iniciais, mais atitude e taxa.
cfg.ci_desafio = [
%   h0     x0    y0    u_B  v_B  vz0     phi°  theta° q°/s
    400,   15,  -10,   0,    0,   -50,   -3,    3,    0.5;
    400,  -12,   15,   0,    0,   -50,    2,   -4,   -0.3;
    500,   20,    0,   2,    0,   -55,    0,    5,    0;
    500,    0,  -18,   0,   -2,   -55,   -5,    0,    0;
    600,   10,   10,   1,   -1,   -50,    3,   -3,    1.0;
    600,  -15,  -12,  -2,    2,   -50,   -4,    4,   -0.5;
    800,    5,   -5,   3,    0,   -50,   -2,    2,    0;
    800,   -8,    8,   0,    3,   -45,    4,   -2,    0.3;
   1000,   18,   15,   1.5,  1.5, -50,   -3,   -3,    0.5;
   1000,  -20,  -15,  -1,   -1,   -50,    5,    5,   -1.0;
    400,    0,    0,   0,    0,   -50,   -1,    1,    0;
    800,   12,   -8,  -2,    2,   -50,    0,    0,    0;
];

cfg.N_CI = size(cfg.ci_desafio, 1);   % número de cenários

% ── Pré-gerar trajetórias ──────────────────────────────────────────────
%   As trajetórias de referência (S-curve) dependem só da CI, não da
%   sintonia. Por isso são geradas UMA vez aqui e reaproveitadas em todas
%   as avaliações.
fprintf('  Pré-gerando trajetórias S-curve (%d CIs)...', cfg.N_CI);
t_traj = tic;

cfg.trajs  = cell(cfg.N_CI, 1);          % uma trajetória por CI
cfg.x0_set = zeros(13, cfg.N_CI);        % vetor de estado inicial (13x1) por CI
traj_falhou = false(cfg.N_CI, 1);        % flag de CIs sem trajetória viável

for i = 1:cfg.N_CI
    ci = cfg.ci_desafio(i, :);
    % Converte ângulos de graus (na tabela) para radianos (no estado).
    phi0_i   = deg2rad(ci(7));
    theta0_i = deg2rad(ci(8));
    q0_i     = deg2rad(ci(9));
    
    % Gera a S-curve de referência para esta CI a partir das posições e
    % velocidades iniciais e da atitude inicial.
    traj_i = gerar_traj_rapido(ci(1), ci(6), ci(2), ci(3), ci(4), ci(5), ...
                                phi0_i, theta0_i, p);
    
    % Algumas CIs podem ser geometricamente inviáveis para a S-curve;
    % nesse caso a CI é marcada e descartada mais abaixo.
    if isempty(traj_i)
        traj_falhou(i) = true;
        fprintf('\n    [AVISO] CI %d: S-curve impossível', i);
        continue;
    end
    
    % Armazena a trajetória e monta o estado inicial completo (13 estados):
    % [x y z | u v w | phi theta psi | p q r | massa].
    cfg.trajs{i} = traj_i;
    cfg.x0_set(:, i) = [ci(2); ci(3); ci(1); ci(4); ci(5); ci(6); ...
                         phi0_i; theta0_i; 0; 0; q0_i; 0; p.m0];
end

% Remove do conjunto as CIs cuja trajetória falhou, mantendo cfg coerente.
if any(traj_falhou)
    boas = ~traj_falhou;
    cfg.ci_desafio = cfg.ci_desafio(boas, :);
    cfg.trajs      = cfg.trajs(boas);
    cfg.x0_set     = cfg.x0_set(:, boas);
    cfg.N_CI       = sum(boas);
end

fprintf(' OK (%.1fs, %d trajetórias)\n', toc(t_traj), cfg.N_CI);


% =========================================================================
% =========================================================================
%  ESTÁGIO 1
%  EXPLORAÇÃO VIA LATIN HYPERCUBE SAMPLING
%
% =========================================================================
% =========================================================================

fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║  ESTÁGIO 1 — EXPLORAÇÃO LHS (%d pontos, %d sims cada)      \n', cfg.N_LHS, cfg.N_CI);
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% ═══════════════════════════════════════════════════════════════════════
%  3. GERAR PONTOS LHS EM LOG-SPACE
% ═══════════════════════════════════════════════════════════════════════
%   Mapeamos pra log-space: knob = lb * (ub/lb)^u,  u ~ U[0,1]
%   Esse mapeamento distribui os pontos uniformemente em ordens de grandeza

rng(42, 'twister');   % reprodutibilidade

% Critério 'maximin' maximiza a distância mínima entre pontos, melhorando
% ainda mais o espalhamento da amostra.
U = lhsdesign(cfg.N_LHS, cfg.N_KNOBS, 'Criterion', 'maximin');

% Mapear [0,1] → [lb, ub] em escala log, knob a knob.
X_lhs = zeros(cfg.N_LHS, cfg.N_KNOBS);
for j = 1:cfg.N_KNOBS
    X_lhs(:, j) = cfg.lb(j) * (cfg.ub(j) / cfg.lb(j)).^U(:, j);
end

% Eco da configuração da amostragem (diagnóstico).
fprintf('  LHS: %d pontos em %dD (log-space, seed=42)\n', cfg.N_LHS, cfg.N_KNOBS);
fprintf('  Bounds:\n');
for j = 1:cfg.N_KNOBS
    fprintf('    %-8s: [%.2f, %.1f]\n', cfg.knob_names{j}, cfg.lb(j), cfg.ub(j));
end

%% ═══════════════════════════════════════════════════════════════════════
%  4. AVALIAR TODOS OS PONTOS LHS
% ═══════════════════════════════════════════════════════════════════════
%   Loop principal do estágio 1: cada um dos N_LHS pontos é convertido em
%   matrizes Q/R, simulado sobre todas as CIs e resumido num custo escalar.
%   É a parte mais cara do script (N_LHS x N_CI simulações).

% Pré-alocação dos vetores de resultado (custo + métricas p95 por ponto).
custos_lhs = NaN(cfg.N_LHS, 1);
metricas_lhs = struct('p95_vz', NaN(cfg.N_LHS,1), ...
                      'p95_lat', NaN(cfg.N_LHS,1), ...
                      'p95_att', NaN(cfg.N_LHS,1), ...
                      'p95_theta', NaN(cfg.N_LHS,1), ...
                      'n_ok', zeros(cfg.N_LHS,1));

fprintf('\n  Avaliando %d pontos LHS...\n', cfg.N_LHS);
t_lhs_start = tic;

for k = 1:cfg.N_LHS
    t_k = tic;
    
    % Montar table row (mesmo formato que bayesopt espera)
    % Usar o mesmo formato no LHS e no BO permite reusar avaliar_candidato.
    params_k = array2table(X_lhs(k,:), 'VariableNames', cfg.knob_names);
    
    % Avalia o candidato: dispara o mini Monte-Carlo e devolve custo+métricas.
    [custos_lhs(k), met_k] = avaliar_candidato(params_k, cfg);
    
    % Guarda as métricas detalhadas deste ponto para a análise posterior.
    metricas_lhs.p95_vz(k)    = met_k.p95_vz;
    metricas_lhs.p95_lat(k)   = met_k.p95_lat;
    metricas_lhs.p95_att(k)   = met_k.p95_att;
    metricas_lhs.p95_theta(k) = met_k.p95_theta;
    metricas_lhs.n_ok(k)      = met_k.n_ok;
    
    % Log periódico de progresso (a cada 10 pontos e no primeiro).
    if mod(k, 10) == 0 || k == 1
        melhor_ate_agr = min(custos_lhs(1:k));
        fprintf('  [%3d/%d] custo=%.1f  melhor=%.1f  (%.0fs/eval)\n', ...
            k, cfg.N_LHS, custos_lhs(k), melhor_ate_agr, toc(t_k));
    end
end

tempo_lhs = toc(t_lhs_start);
fprintf('\n  Estágio 1 concluído: %.1f min (%d avaliações)\n', tempo_lhs/60, cfg.N_LHS);

%% ═══════════════════════════════════════════════════════════════════════
%  5. ANÁLISE DE SENSIBILIDADE
% ═══════════════════════════════════════════════════════════════════════
%   A partir dos 200 pontos avaliados, decide-se quais knobs importam e
%   quais já convergiram. Essa análise é o que conecta o estágio 1 ao 2:
%   ela define quem fixar e quão estreitos serão os novos bounds.

fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║  ANÁLISE DE SENSIBILIDADE                                    ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

% ── 5a. Correlação de Spearman ─────────────────────────────────────────
%   Mede a influência monotônica de cada knob no custo.
%   |rho| alto → knob importa; ~0 → knob irrelevante.
%   Spearman (e não Pearson) porque a relação knob→custo é não-linear, mas
%   tipicamente monotônica em escala log.

% Considera só os pontos com custo válido (descarta eventuais NaN).
validos = ~isnan(custos_lhs);
X_val   = X_lhs(validos, :);
C_val   = custos_lhs(validos);

rho_spearman = zeros(cfg.N_KNOBS, 1);   % coeficiente de correlação
p_valor      = zeros(cfg.N_KNOBS, 1);   % significância estatística

% Correlaciona cada knob individualmente contra o custo.
for j = 1:cfg.N_KNOBS
    [rho_spearman(j), p_valor(j)] = corr(X_val(:,j), C_val, 'Type', 'Spearman');
end

% Tabela de influência no Command Window.
fprintf('  CORRELAÇÃO DE SPEARMAN (knob vs custo):\n');
fprintf('  %-8s  |  rho      |  |rho|    |  p-valor  |  Influência\n', 'Knob');
fprintf('  ---------+----------+----------+-----------+-------------\n');

for j = 1:cfg.N_KNOBS
    abs_rho = abs(rho_spearman(j));
    % Classificação qualitativa por faixas de |rho| (apenas para leitura).
    if abs_rho > 0.3
        nivel = 'FORTE';
    elseif abs_rho > 0.15
        nivel = 'moderada';
    else
        nivel = 'fraca';
    end
    fprintf('  %-8s  | %+.3f    |  %.3f   |  %.4f   |  %s\n', ...
        cfg.knob_names{j}, rho_spearman(j), abs_rho, p_valor(j), nivel);
end

% ── 5b. Análise de cluster dos TOP-N ──────────────────────────────────
%   Pega os N_TOP melhores pontos e calcula estatísticas por knob.
%   CV baixo → knob convergiu → pode fixar.
%   A correlação (5a) diz se o knob importa; o cluster (5b) diz se já
%   sabemos QUAL valor ele deve ter. Só fixamos quem já convergiu.

% Ordena por custo crescente e seleciona o cluster dos melhores pontos.
[~, idx_sort] = sort(custos_lhs);
idx_top = idx_sort(1:min(cfg.N_TOP, sum(validos)));
X_top   = X_lhs(idx_top, :);
C_top   = custos_lhs(idx_top);

fprintf('\n  CLUSTER DOS TOP-%d MELHORES (custo %.1f a %.1f):\n', ...
    length(idx_top), C_top(1), C_top(end));
fprintf('  %-8s  |  mediana   |  p5        |  p95       |  CV     |  Status\n', 'Knob');
fprintf('  ---------+-----------+------------+------------+---------+-----------\n');

% Estatísticas por knob, computadas sobre o cluster dos melhores.
knob_mediana = zeros(cfg.N_KNOBS, 1);   % valor central (em linear)
knob_cv      = zeros(cfg.N_KNOBS, 1);   % coef. de variação (em log)
knob_p5      = zeros(cfg.N_KNOBS, 1);   % percentil 5 (em linear)
knob_p95     = zeros(cfg.N_KNOBS, 1);   % percentil 95 (em linear)
knob_fixo    = false(cfg.N_KNOBS, 1);   % decisão fixar/otimizar

for j = 1:cfg.N_KNOBS
    vals = X_top(:, j);
    
    % Estatísticas em log-space
    % Como os knobs vivem em escala log, o CV é calculado sobre log10 dos
    % valores; isso evita que knobs com média grande pareçam artificialmente
    % "convergidos" só por causa da escala.
    log_vals = log10(vals);
    mu_log   = mean(log_vals);
    std_log  = std(log_vals);
    cv_log   = std_log / abs(mu_log);   % CV em log-space
    
    knob_mediana(j) = median(vals);
    knob_p5(j)      = prctile(vals, 5);
    knob_p95(j)     = prctile(vals, 95);
    knob_cv(j)      = cv_log;
    % Dispersão baixa entre os melhores => knob convergiu => fixar.
    knob_fixo(j)    = cv_log < cfg.CV_THRESHOLD;
    
    if knob_fixo(j)
        status = 'FIXAR';
    else
        status = 'otimizar';
    end
    
    fprintf('  %-8s  | %8.3f  | %8.3f   | %8.3f   |  %.3f  |  %s\n', ...
        cfg.knob_names{j}, knob_mediana(j), knob_p5(j), knob_p95(j), ...
        knob_cv(j), status);
end

% Contagem de quantos knobs vão para o estágio 2.
n_fixos = sum(knob_fixo);
n_livres = cfg.N_KNOBS - n_fixos;

fprintf('\n  Resultado: %d knobs fixados, %d livres para estágio 2\n', n_fixos, n_livres);

% ── 5c. Gráficos ──────────────────────────────────────────────────────
%   Painel de diagnóstico: barras de influência + dispersões knob-vs-custo.
%   Serve para inspeção visual e como figura de apoio no relatório do TCC.

figure('Name', 'Stage 1 — Sensibility Analysis', ...
       'NumberTitle', 'off', 'Position', [50 50 1400 800]);

% Barplot de |rho|
% Quão fortemente cada knob influencia o custo, com limiares de referência.
subplot(2, 4, 1);
bar(abs(rho_spearman), 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', cfg.knob_names, 'XTickLabelRotation', 45);
ylabel('|Spearman \rho|');
title('Cost Influence');
yline(0.15, 'r--', 'Moderate threshold');
yline(0.30, 'r-',  'Strong threshold');
grid on;

% Scatter de cada knob vs custo (top 100 pontos pra não poluir)
% Cada subgráfico mostra a "dependência parcial" aproximada do custo em
% relação a um knob; a cor reforça o nível de custo de cada ponto.
idx_plot = idx_sort(1:min(100, length(idx_sort)));
for j = 1:cfg.N_KNOBS
    subplot(2, 4, 1 + j);
    scatter(X_lhs(idx_plot, j), custos_lhs(idx_plot), 15, ...
            custos_lhs(idx_plot), 'filled', 'MarkerFaceAlpha', 0.6);
    set(gca, 'XScale', 'log');   % eixo log: coerente com a amostragem
    xlabel(cfg.knob_names{j});
    ylabel('Cost');
    title(sprintf('%s (\\rho=%.2f)', cfg.knob_names{j}, rho_spearman(j)));
    colormap(gca, flipud(hot));
    
    % Marcar a mediana do TOP-N
    % Linha azul: valor central dos melhores (candidato a valor fixo).
    xline(knob_mediana(j), 'b-', 'LineWidth', 1.5);
    
    % Marcar bounds do estágio 2 se não fixado
    % Linhas verdes: faixa [p5, p95] que dará origem aos bounds refinados.
    if ~knob_fixo(j)
        xline(knob_p5(j),  'g--', 'p5');
        xline(knob_p95(j), 'g--', 'p95');
    end
    grid on;
end

sgtitle('Estágio 1 — LHS: Sensibilidade dos 7 Knobs');
drawnow;

%% ═══════════════════════════════════════════════════════════════════════
%  6. PREPARAR ESTÁGIO 2
% ═══════════════════════════════════════════════════════════════════════
%   Traduz as decisões da análise de sensibilidade em parâmetros concretos
%   do BO: valores congelados para os knobs fixos e bounds estreitos para
%   os knobs livres.

% ── Valores fixados (mediana dos TOP-N) ────────────────────────────────
%   Knobs congelados assumem a mediana do cluster dos melhores — o valor
%   mais representativo das boas soluções.
knob_fixo_val = knob_mediana;   % valor em que fixar

% ── Bounds apertados para knobs livres ─────────────────────────────────
%   Usa [p5, p95] dos TOP-N com margem de 2× pra cada lado (log-space)
%   A margem evita "amarrar" o BO exatamente à amostra do LHS, deixando
%   espaço para ele explorar um pouco além do que o LHS já viu.
margem_bound = 2.0;   % fator multiplicativo

% Parte dos bounds originais e só estreita os knobs que permanecem livres.
lb_refinado = cfg.lb;
ub_refinado = cfg.ub;

for j = 1:cfg.N_KNOBS
    if ~knob_fixo(j)
        % Aperta os limites em torno da faixa boa, sem estourar os bounds
        % globais originais (max/min garantem que nunca se expanda além).
        lb_refinado(j) = max(cfg.lb(j), knob_p5(j)  / margem_bound);
        ub_refinado(j) = min(cfg.ub(j), knob_p95(j) * margem_bound);
        
        % Garantir que o range não ficou degenerado
        % Se a faixa ficou estreita demais (< meia década), recentra em
        % torno da mediana com um fator 3× para o BO ter espaço de manobra.
        if ub_refinado(j) / lb_refinado(j) < 2
            centro = knob_mediana(j);
            lb_refinado(j) = max(cfg.lb(j), centro / 3);
            ub_refinado(j) = min(cfg.ub(j), centro * 3);
        end
    end
end

% Resumo da configuração do estágio 2 (fixos vs. livres e seus bounds).
fprintf('\n  CONFIGURAÇÃO DO ESTÁGIO 2:\n');
fprintf('  %-8s  |  Status      |  Valor/Bounds\n', 'Knob');
fprintf('  ---------+--------------+-----------------------------\n');
for j = 1:cfg.N_KNOBS
    if knob_fixo(j)
        fprintf('  %-8s  |  FIXO        |  %.4f (mediana TOP-%d)\n', ...
            cfg.knob_names{j}, knob_fixo_val(j), cfg.N_TOP);
    else
        fprintf('  %-8s  |  LIVRE       |  [%.3f, %.3f] (era [%.2f, %.1f])\n', ...
            cfg.knob_names{j}, lb_refinado(j), ub_refinado(j), cfg.lb(j), cfg.ub(j));
    end
end

fprintf('\n  Dimensionalidade: %dD → %dD\n', cfg.N_KNOBS, n_livres);

% ── Salvar resultados do estágio 1 ────────────────────────────────────
%   Persistir tudo permite reanalisar o estágio 1 sem ter que reexecutá-lo
%   (que leva dezenas de minutos).
save('resultados_estagio1.mat', ...
    'X_lhs', 'custos_lhs', 'metricas_lhs', ...
    'rho_spearman', 'p_valor', ...
    'X_top', 'C_top', 'idx_top', ...
    'knob_mediana', 'knob_cv', 'knob_p5', 'knob_p95', ...
    'knob_fixo', 'knob_fixo_val', ...
    'lb_refinado', 'ub_refinado', 'cfg');

fprintf('\n  Estágio 1 salvo em resultados_estagio1.mat\n');

% =========================================================================
% =========================================================================
%  ESTÁGIO 2
%  REFINAMENTO VIA BAYESIAN OPTIMIZATION
%
% =========================================================================
% =========================================================================

fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║  ESTÁGIO 2 — REFINAMENTO BAYESOPT (%dD, %d evals)          \n', n_livres, cfg.N_BAYESOPT);
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

% Caso de borda: se a análise fixou TODOS os knobs, não há o que otimizar e
% a solução já é a combinação de medianas dos melhores pontos.
if n_livres == 0
    fprintf('  Todos os knobs convergiram no estágio 1!\n');
    fprintf('  Usando medianas do TOP-%d diretamente.\n', cfg.N_TOP);
    best_final = knob_fixo_val;
else
    % ── Montar variáveis otimizáveis (só as livres) ────────────────────
    %   Cada knob livre vira uma optimizableVariable com bounds refinados e
    %   transformação log, coerente com a escala usada o tempo todo.
    vars_livres = [];
    idx_livres  = find(~knob_fixo);
    
    for jj = 1:length(idx_livres)
        j = idx_livres(jj);
        vars_livres = [vars_livres; ...
            optimizableVariable(cfg.knob_names{j}, ...
                [lb_refinado(j), ub_refinado(j)], 'Transform', 'log')]; %#ok<AGROW>
    end
    
    fprintf('  Variáveis livres: ');
    fprintf('%s ', cfg.knob_names{idx_livres});
    fprintf('\n');
    
    % ── Chute inicial: melhor ponto do LHS ─────────────────────────────
    %   Aproveita o melhor candidato já encontrado como ponto de partida do
    %   BO (InitialX), acelerando a convergência.
    melhor_lhs = X_lhs(idx_sort(1), :);
    chute_table = array2table(melhor_lhs(idx_livres), ...
                              'VariableNames', cfg.knob_names(idx_livres));
    
    fprintf('  Chute (melhor LHS): ');
    for jj = 1:length(idx_livres)
        fprintf('%s=%.3f ', cfg.knob_names{idx_livres(jj)}, melhor_lhs(idx_livres(jj)));
    end
    fprintf('\n');
    
    % ── Closure: knobs fixos embutidos ─────────────────────────────────
    %   O BO só enxerga os knobs livres; este wrapper recombina os fixos
    %   antes de cada avaliação, mantendo a função-objetivo consistente.
    fun_bo = @(params_livre) avaliar_com_fixos(params_livre, ...
        knob_fixo, knob_fixo_val, cfg);
    
    % ── Rodar bayesopt ─────────────────────────────────────────────────
    %   Otimização bayesiana com 'expected-improvement-plus'. IsObjectiveDeterministic=true porque o custo
    %   é reprodutível (CIs e seed fixas). ExplorationRatio=0.5 dá um peso
    %   maior à exploração, evitando travar num ótimo local cedo demais.
    fprintf('\n  Iniciando Bayesian Optimization...\n\n');
    t_bo_start = tic;
    
    resultados_bo = bayesopt(fun_bo, vars_livres, ...
        'MaxObjectiveEvaluations', cfg.N_BAYESOPT, ...
        'NumSeedPoints',           cfg.N_SEED_BO, ...
        'InitialX',                chute_table, ...
        'IsObjectiveDeterministic', true, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ExplorationRatio',        0.5, ...
        'UseParallel',             false, ...
        'PlotFcn',                 {@plotObjectiveModel, @plotMinObjective}, ...
        'Verbose',                 1);
    
    tempo_bo = toc(t_bo_start);
    
    % ── Reconstruir vetor completo de 7 knobs ──────────────────────────
    %   Combina o melhor ponto do BO (knobs livres) com os valores fixos
    %   para formar a sintonia final de 7 knobs.
    best_bo = resultados_bo.XAtMinObjective;
    best_final = knob_fixo_val;   % começa com os fixos
    
    for jj = 1:length(idx_livres)
        j = idx_livres(jj);
        best_final(j) = best_bo.(cfg.knob_names{j});
    end
end

%% ═══════════════════════════════════════════════════════════════════════
%  7. RESULTADO FINAL
% ═══════════════════════════════════════════════════════════════════════
%   Constrói as matrizes finais, reavalia a sintonia escolhida para gerar o
%   relatório de desempenho e imprime um bloco pronto para colar no init_mpc.

% Montar Q e R finais a partir dos 7 knobs ótimos.
params_final = array2table(best_final', 'VariableNames', cfg.knob_names);
[Q_opt, R_opt] = montar_QR(params_final, cfg);

% Disponibiliza Q_opt/R_opt no workspace base para uso posterior.
assignin('base', 'Q_mpc', Q_opt);
assignin('base', 'R_mpc', R_opt);

% Re-avaliar pra relatório
% Roda o mini Monte-Carlo uma última vez na solução final para obter as
% métricas que serão exibidas e comparadas com os requisitos.
[custo_final, met_final] = avaliar_candidato(params_final, cfg);

% Painel de resultado final.
fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║              RESULTADO FINAL (2 ESTÁGIOS)                    ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════╣\n');
fprintf('║  Score final: %.2f                                          \n', custo_final);
fprintf('╠═══════════════════════════════════════════════════════════════╣\n');
fprintf('║  MULTIPLICADORES:                                            ║\n');
% Lista cada knob final, indicando se veio fixado (E1) ou refinado (E2).
for j = 1:cfg.N_KNOBS
    if knob_fixo(j)
        tag = '(fixado E1)';
    else
        tag = '(refinado E2)';
    end
    fprintf('║    %-8s = %8.4f   %s\n', cfg.knob_names{j}, best_final(j), tag);
end
fprintf('╠═══════════════════════════════════════════════════════════════╣\n');

% Confronta cada métrica p95 com seu requisito, marcando OK/FALHA.
if met_final.n_ok > 0
    check = @(v, r) sel_str(v <= r, 'OK', 'FALHA');
    fprintf('║  MÉTRICAS p95 (teste inicial e raso só para expectativa):  ║\n');
    fprintf('║    Erro lateral:  %5.2f m   (req %.1f)  %s                  \n', met_final.p95_lat, cfg.REQ_LAT, check(met_final.p95_lat, cfg.REQ_LAT));
    fprintf('║    |vz| td:      %5.2f m/s (req %.1f)  %s                  \n', met_final.p95_vz, cfg.REQ_VZ, check(met_final.p95_vz, cfg.REQ_VZ));
    fprintf('║    |att| td:     %5.2f °   (req %.1f)  %s                  \n', met_final.p95_att, cfg.REQ_ATT, check(met_final.p95_att, cfg.REQ_ATT));
    fprintf('║    theta_max:    %5.2f °   (req %.1f)  %s                  \n', met_final.p95_theta, cfg.REQ_THETA_MAX, check(met_final.p95_theta, cfg.REQ_THETA_MAX));
    fprintf('║    Pousos:       %d/%d                                      \n', met_final.n_ok, cfg.N_CI);
end

fprintf('╚═══════════════════════════════════════════════════════════════╝\n');

% ── Bloco pra colar no init_mpc ────────────────────────────────────────
%   Emite o código-fonte dos knobs e das matrizes Q/R já com os comentários
%   de mapeamento estado-a-estado, pronto para ser transplantado no setup.
fprintf('\n  ── COLAR NO init_mpc.m ────────────────────────────────\n');
for j = 1:cfg.N_KNOBS
    fprintf('  %-8s = %.4f;\n', cfg.knob_names{j}, best_final(j));
end
fprintf('\n');
fprintf('  Q_mpc = diag([...\n');
fprintf('      k_pos  / (5.0)^2,           ... %%  1. x_I\n');
fprintf('      k_pos  / (5.0)^2,           ... %%  2. y_I\n');
fprintf('      k_h    / (2.0)^2,           ... %%  3. z_I\n');
fprintf('      k_vel  / (3.0)^2,           ... %%  4. u_B\n');
fprintf('      k_vel  / (3.0)^2,           ... %%  5. v_B\n');
fprintf('      k_vz   / (1.0)^2,           ... %%  6. w_B\n');
fprintf('      k_ang  / (deg2rad(5))^2,    ... %%  7. phi\n');
fprintf('      k_ang  / (deg2rad(5))^2,    ... %%  8. theta\n');
fprintf('      k_rate / (deg2rad(5))^2,    ... %%  9. p_B\n');
fprintf('      k_rate / (deg2rad(5))^2     ... %% 10. q_B\n');
fprintf('  ]);\n');
fprintf('  R_mpc = diag([...\n');
fprintf('      1/(p.T_max - p.T_min)^2,    ... %% T\n');
fprintf('      k_tvc / (p.delta_max)^2,     ... %% delta_y\n');
fprintf('      k_tvc / (p.delta_max)^2      ... %% delta_z\n');
fprintf('  ]);\n');
fprintf('  ──────────────────────────────────────────────────────\n');

%% ═══════════════════════════════════════════════════════════════════════
%  8. SALVAR TUDO
% ═══════════════════════════════════════════════════════════════════════
%   Consolida em um único .mat o resultado final.

save('resultados_2estagios.mat', ...
    'best_final', 'custo_final', 'met_final', ...
    'Q_opt', 'R_opt', ...
    'X_lhs', 'custos_lhs', 'metricas_lhs', ...
    'rho_spearman', 'knob_fixo', 'knob_fixo_val', ...
    'knob_mediana', 'knob_cv', 'knob_p5', 'knob_p95', ...
    'cfg');

% Os resultados do BO só existem se houve knobs livres; anexa se for o caso.
if n_livres > 0
    save('resultados_2estagios.mat', 'resultados_bo', '-append');
end

fprintf('\n  Tudo salvo em resultados_2estagios.mat\n');
fprintf('  Próximo passo: validar com Monte Carlo completo.\n\n');


%% ═════════════════════════════════════════════════════════════════════
%  FUNÇÕES LOCAIS
% ═════════════════════════════════════════════════════════════════════
%   Funções auxiliares usadas tanto pelo LHS quanto pelo BO. Concentrar a
%   lógica de avaliação aqui garante que os dois estágios meçam o custo
%   exatamente da mesma forma.

function [custo, met] = avaliar_candidato(params, cfg)
% AVALIAR_CANDIDATO  Roda mini-MC com 12 CIs e retorna custo + métricas.
%
%   É o coração da função-objetivo: converte os knobs em Q/R, simula todas
%   as CIs (em paralelo ou serial) e agrega o desempenho num único escalar.
%
%   params — table row com 7 knobs (k_pos, k_vel, ..., k_tvc)
%   cfg    — struct de configuração
%   custo  — escalar (menor = melhor)
%   met    — struct com p95_vz, p95_lat, p95_att, p95_theta, n_ok

    % Constrói as matrizes de ponderação a partir dos knobs.
    [Q_test, R_test] = montar_QR(params, cfg);
    
    % Publica Q/R no workspace base (o modelo Simulink lê dali).
    assignin('base', 'Q_mpc', Q_test);
    assignin('base', 'R_mpc', R_test);
    
    % Vetores para coletar as métricas de cada cenário.
    N_CI = cfg.N_CI;
    err_lat   = NaN(N_CI, 1);   % erro lateral no toque [m]
    vz_td     = NaN(N_CI, 1);   % |velocidade vertical| no toque [m/s]
    att_td    = NaN(N_CI, 1);   % erro de atitude no toque [graus]
    theta_max = NaN(N_CI, 1);   % pico de theta na descida [graus]
    pousou    = false(N_CI, 1); % flag de pouso bem-sucedido
    
    if cfg.USAR_PARSIM
        % ── Caminho paralelo: monta um SimulationInput por CI ──────────
        % Cada entrada carrega suas próprias Q/R, estado inicial, trajetória
        % e tempo de parada (duração da trajetória + 5 s de folga).
        simIn(1:N_CI) = Simulink.SimulationInput(cfg.modelo);
        for i = 1:N_CI
            t_max_i = cfg.trajs{i}.T + 5;
            simIn(i) = simIn(i).setVariable('Q_mpc', Q_test);
            simIn(i) = simIn(i).setVariable('R_mpc', R_test);
            simIn(i) = simIn(i).setVariable('x0_sim', cfg.x0_set(:,i));
            simIn(i) = simIn(i).setVariable('traj', cfg.trajs{i});
            simIn(i) = simIn(i).setModelParameter('StopTime', num2str(t_max_i));
            % PreSimFcn força as variáveis no workspace de cada worker.
            simIn(i) = simIn(i).setPreSimFcn(@forcar_workspace);
        end
        
        % Dispara todas as simulações em paralelo.
        simOut = parsim(simIn, ...
            'ShowProgress', 'off', ...
            'TransferBaseWorkspaceVariables', 'on');
        
        % Extrai as métricas de cada saída.
        for i = 1:N_CI
            [err_lat(i), vz_td(i), att_td(i), theta_max(i), pousou(i)] = ...
                extrair_metricas(simOut(i));
        end
    else
        % ── Caminho serial: simula uma CI de cada vez ──────────────────
        % Faz backup das variáveis do workspace base para restaurá-las ao
        % final, evitando efeitos colaterais fora desta função.
        x0_backup   = evalin('base', 'x0_sim');
        traj_backup = evalin('base', 'traj');
        
        for i = 1:N_CI
            assignin('base', 'x0_sim', cfg.x0_set(:,i));
            assignin('base', 'traj',   cfg.trajs{i});
            t_max_i = cfg.trajs{i}.T + 5;
            try
                simOut_i = sim(cfg.modelo, 'StopTime', num2str(t_max_i), ...
                              'SrcWorkspace', 'base');
                [err_lat(i), vz_td(i), att_td(i), theta_max(i), pousou(i)] = ...
                    extrair_metricas(simOut_i);
            catch
                % Falha de simulação
            end
        end
        
        % Restaura o workspace base ao estado original.
        assignin('base', 'x0_sim', x0_backup);
        assignin('base', 'traj',   traj_backup);
    end
    
    % Contagem de sucessos e falhas de pouso.
    n_ok    = sum(pousou);
    n_crash = N_CI - n_ok;
    
    % Inicializa a struct de métricas (NaN até haver dados suficientes).
    met.n_ok      = n_ok;
    met.p95_vz    = NaN;
    met.p95_lat   = NaN;
    met.p95_att   = NaN;
    met.p95_theta = NaN;
    
    if n_ok >= 3
        % Com pelo menos 3 pousos, calcula o p95 das métricas sobre os
        % casos bem-sucedidos. O p95 captura o "pior caso típico", que é o
        % que interessa para robustez (mais informativo que a média).
        met.p95_vz    = prctile(vz_td(pousou),     95);
        met.p95_lat   = prctile(err_lat(pousou),    95);
        met.p95_att   = prctile(att_td(pousou),     95);
        met.p95_theta = prctile(theta_max(pousou),  95);
        med_vz        = median(vz_td(pousou));
        
        % Custo ponderado: prioriza velocidade vertical e atitude no toque
        % (mais críticas para um pouso seguro), com pesos menores para erro
        % lateral, theta e a mediana de vz, e penalidade forte por crash.
        custo = 10 * met.p95_vz  ...
              +  3 * met.p95_lat ...
              +  4 * met.p95_att ...
              +  0.5 * met.p95_theta ...
              +  0.5 * med_vz ...
              + 50 * (n_crash / N_CI);
    else
        % Poucos pousos => sintonia ruim. Atribui custo-base alto, ainda
        % proporcional à taxa de crash, para que o otimizador consiga
        % distinguir entre soluções igualmente ruins e fugir delas.
        custo = 500 + 100 * (n_crash / N_CI);
    end
end


function custo = avaliar_com_fixos(params_livre, knob_fixo, knob_fixo_val, cfg)
% AVALIAR_COM_FIXOS  Wrapper que reconstrói os 7 knobs (fixos + livres)
%   e chama avaliar_candidato.
%
%   É a função-objetivo vista pelo bayesopt: ele só fornece os knobs livres,
%   e este wrapper os combina com os fixos antes de avaliar. Também imprime
%   o progresso de cada avaliação do BO.

    % Contador persistente para numerar as avaliações do BO no log.
    persistent eval_count_bo
    if isempty(eval_count_bo); eval_count_bo = 0; end
    eval_count_bo = eval_count_bo + 1;
    t_eval = tic;
    
    % Reconstruir vetor completo
    % Começa com todos os knobs nos valores fixos e sobrescreve apenas os
    % livres com o que o BO propôs nesta iteração.
    vals = knob_fixo_val;   % começa com tudo fixo
    
    nomes_livres = params_livre.Properties.VariableNames;
    for jj = 1:length(nomes_livres)
        idx = find(strcmp(cfg.knob_names, nomes_livres{jj}));
        vals(idx) = params_livre.(nomes_livres{jj});
    end
    
    % Avalia o vetor completo de 7 knobs reusando avaliar_candidato.
    params_full = array2table(vals', 'VariableNames', cfg.knob_names);
    [custo, met] = avaliar_candidato(params_full, cfg);
    
    % Log da avaliação: formato detalhado se pousou, compacto se "crashou".
    t_el = toc(t_eval);
    if met.n_ok >= 3
        fprintf('  [BO %3d] custo=%.1f | ok=%d/%d | p95: vz=%.2f lat=%.2f att=%.1f theta=%.1f | %.0fs\n', ...
            eval_count_bo, custo, met.n_ok, cfg.N_CI, ...
            met.p95_vz, met.p95_lat, met.p95_att, met.p95_theta, t_el);
    else
        fprintf('  [BO %3d] custo=%.1f | ok=%d/%d (crash) | %.0fs\n', ...
            eval_count_bo, custo, met.n_ok, cfg.N_CI, t_el);
    end
end


function [Q, R] = montar_QR(params, cfg)
% MONTAR_QR  Constrói Q_mpc e R_mpc a partir dos 7 knobs.
%   Aplica a regra de Bryson: cada peso = knob / (valor_máximo)^2. Knobs que
%   se repetem (ex.: k_pos em x e y) refletem a simetria do problema; assim,
%   7 knobs parametrizam os 10 termos de Q e os 3 de R.
    Q = diag([
        params.k_pos  / cfg.bryson.pos^2     % x_I  — posição lateral X
        params.k_pos  / cfg.bryson.pos^2     % y_I  — posição lateral Y
        params.k_h    / cfg.bryson.h^2       % z_I  — altitude
        params.k_vel  / cfg.bryson.vel^2     % u_B  — velocidade X (corpo)
        params.k_vel  / cfg.bryson.vel^2     % v_B  — velocidade Y (corpo)
        params.k_vz   / cfg.bryson.vz^2      % w_B  — velocidade vertical
        params.k_ang  / cfg.bryson.ang^2     % phi  — rolagem
        params.k_ang  / cfg.bryson.ang^2     % theta— arfagem
        params.k_rate / cfg.bryson.rate^2    % p_B  — taxa de rolagem
        params.k_rate / cfg.bryson.rate^2    % q_B  — taxa de arfagem
    ]);
    R = diag([
        cfg.R_base_T                 % T       — empuxo (peso base fixo)
        params.k_tvc * cfg.R_base_tvc % delta_y — deflexão TVC eixo y
        params.k_tvc * cfg.R_base_tvc % delta_z — deflexão TVC eixo z
    ]);
end


function simIn = forcar_workspace(simIn)
% FORCAR_WORKSPACE  PreSimFcn para parsim.
%   Executada por cada worker logo antes da simulação: copia as variáveis
%   do SimulationInput para o workspace base do worker, garantindo que o
%   modelo encontre Q_mpc, R_mpc, x0_sim e traj corretos.
    assignin('base', 'Q_mpc',  simIn.getVariable('Q_mpc'));
    assignin('base', 'R_mpc',  simIn.getVariable('R_mpc'));
    assignin('base', 'x0_sim', simIn.getVariable('x0_sim'));
    assignin('base', 'traj',   simIn.getVariable('traj'));
end


function [err_lat, vz_td, att_td, theta_max, pousou] = extrair_metricas(simOut)
% EXTRAIR_METRICAS  Extrai métricas de pouso de um SimulationOutput.
%   Lê o histórico de estados da simulação e calcula as métricas de toque.
%   É robusta a falhas: qualquer erro devolve NaN/"não pousou" em vez de
%   interromper a otimização.

    % Defaults para o caso de a simulação ter falhado ou ser inválida.
    err_lat   = NaN;
    vz_td     = NaN;
    att_td    = NaN;
    theta_max = NaN;
    pousou    = false;
    
    % Se o Simulink reportou erro, aborta a extração silenciosamente.
    if isprop(simOut, 'ErrorMessage') && ~isempty(simOut.ErrorMessage)
        return;
    end
    
    try
        % Histórico de estados como matriz [n_amostras x 13].
        x_data = squeeze(simOut.x_sim.Data)';
        if size(x_data, 1) < 10;  return;  end   % poucos pontos => inválido
        
        x_f = x_data(end, :);   % estado no instante do toque
        
        % Métricas de toque a partir do estado final.
        err_lat   = sqrt(x_f(1)^2 + x_f(2)^2);            % distância lateral
        vz_td     = abs(x_f(6));                          % |vz| no toque
        att_td    = rad2deg(sqrt(x_f(7)^2 + x_f(8)^2));   % erro de atitude
        theta_max = rad2deg(max(abs(x_data(:, 8))));      % pico de theta
        
        % Critério de pouso bem-sucedido: próximo do solo, descendo devagar
        % e com atitude controlada.
        pousou = abs(x_f(3)) < 10 && vz_td < 10 && att_td < 30;
    catch
        % Qualquer falha na extração mantém os defaults (NaN / não pousou).
    end
end


function r = sel_str(cond, v_true, v_false)
% SEL_STR  Seleção de string (ternário).
%   Pequeno utilitário que emula o operador ternário inexistente em MATLAB,
%   usado para imprimir "OK"/"FALHA" no relatório final.
    if cond; r = v_true; else; r = v_false; end
end