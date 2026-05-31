% =========================================================================
% INIT_EKF.M  —  Tuning do Filtro de Kalman Estendido (EKF)
%%
%   1. Estado expandido de 13 → 16 (acrescenta b_gx, b_gy, b_gz).
%      Resolve a inconsistência arquitetural em que o Sensor_Model injeta
%      bias constante nos giroscópios mas o EKF tratava como ruído branco.
%      Com o bias como estado, o filtro o estima online e remove o "drift"
%      sistemático da medida de taxa angular.
%
%   2. Massa SAI do update: y_meas(13) (massa) NÃO é mais usada na correção.
%      A massa é apenas integrada pela equação de Tsiolkovsky no Predict.
%      Modelagem fisicamente correta — foguete real não mede massa em voo.
%
%   3. Histerese na transição de regime (radar lock):
%      - Entra em modo pouso quando h < 150m  (h_radar_lock)
%      - Sai do modo pouso quando h > 170m   (h_radar_unlock)
%      Evita chattering do filtro caso a estimativa oscile em torno de 150m.
%
%   4. O EKF usa Joseph form para atualizar P (mais robusto numericamente).
%
%   ── INTERFACE COM O SIMULINK ─────────────────────────────────────────
%
%   Tudo permanece compatível com a fiação existente do .slx:
%     - Entrada do bloco Kalman_Filter: y_meas (13×1) 
%     - Saída do bloco Kalman_Filter:   x_hat  (13×1)
%     - Parâmetros lidos do workspace:  P0_ekf, Q_ekf, R_ekf_voo,
%                                       R_ekf_pouso, h_radar_lock
%   O bloco internamente passa de 13 para 16 estados, mas a saída final
%   continua sendo o vetor [pos; vel; euler; omega; m] que o MPC consome.
%
% Execute antes de rodar a simulação (após init_sim e init_mpc).
% =========================================================================
fprintf('\n=== Inicializando Parâmetros do EKF ===\n');

% ─────────────────────────────────────────────────────────────────────────
% Layout do vetor de estado interno do EKF (16 elementos):
%
%   x_est( 1: 3) = posição inercial [m]
%   x_est( 4: 6) = velocidade no corpo [m/s]
%   x_est( 7: 9) = ângulos de Euler [rad]
%   x_est(10:12) = taxas angulares no corpo [rad/s]
%   x_est(13)    = massa [kg]
%   x_est(14:16) = bias do giroscópio [rad/s]
%
% Layout das medidas (Sensor_Model continua entregando 13×1):
%   y_meas(1:12) — pos, vel, euler, omega         (usadas pelo EKF)
%   y_meas(13)   — massa                          (NÃO usada, só predição)
% ─────────────────────────────────────────────────────────────────────────

% ─────────────────────────────────────────────────────────────────────────
% 1. Covariância Inicial P0 (16×16)
%
%    Para o bias, σ_inicial = 2°/s. Esse valor é GRANDE de propósito:
%    o bias real injetado pelo Sensor_Model é ~0.02 rad/s ≈ 1.15°/s,
%    portanto precisamos que o 1σ inicial seja maior que ele para o filtro
%    aceitar correções significativas no início. Conforme o bias converge,
%    P(14:16,14:16) cai naturalmente.
% ─────────────────────────────────────────────────────────────────────────
sigma_bias_init = deg2rad(2.0);   % 1σ inicial do bias [rad/s] (~ 2°/s)

% Diagonal de P0 = variâncias iniciais (quadrado dos 1σ) de cada estado.
P0_ekf = diag([ ...
    10, 10, 10,                                    ... %  1-3   posição    [m²]
    1,  1,  1,                                     ... %  4-6   velocidade [(m/s)²]
    0.1, 0.1, 0.1,                                 ... %  7-9   atitude    [rad²]
    0.01, 0.01, 0.01,                              ... % 10-12  taxa       [(rad/s)²]
    100,                                           ... % 13     massa      [kg²]
    sigma_bias_init^2, sigma_bias_init^2, sigma_bias_init^2 ... % 14-16 bias
]);

% ─────────────────────────────────────────────────────────────────────────
% 2. Ruído de Processo Q (16×16)
% ─────────────────────────────────────────────────────────────────────────
% Desvios-padrão por grupo de estado (1σ); a variância entra em Q via .^2.
std_model_pos   = 1e-2;             % [m]      — confiança alta na cinemática translacional
std_model_vel   = 1e-2;             % [m/s]    — mesma coisa
std_model_ang   = deg2rad(0.05);    % [rad]    — relaxado vs. v1 (era 0.02°)
std_model_taxa  = deg2rad(0.02);    % [rad/s]  — relaxado vs. v1 (era 0.01°)
std_model_m     = 1.0;              % [kg]     — Tsiolkovsky é exata, mas dt finito
std_model_bias  = deg2rad(0.001);   % [rad/s/(passo^1/2)] — random walk do bias

% Replicação por eixo (×3) e quadratura para formar a diagonal de Q.
Q_ekf = diag([ ...
    std_model_pos   * ones(1,3), ...
    std_model_vel   * ones(1,3), ...
    std_model_ang   * ones(1,3), ...
    std_model_taxa  * ones(1,3), ...
    std_model_m,                 ...
    std_model_bias  * ones(1,3)  ...
].^2);

% ─────────────────────────────────────────────────────────────────────────
% 3. Ruído de Medição R — ALTA ALTITUDE (GPS comum + Barômetro)
%
%    R é mantida com 13 linhas/colunas porque o Sensor_Model continua
%    produzindo y_meas com 13 elementos (compatibilidade do .slx).
%    O EKF internamente extrai R(1:12, 1:12) para o update.
%    A linha 13 (massa) fica como placeholder e é ignorada pelo filtro.
% ─────────────────────────────────────────────────────────────────────────
% 1σ de cada sensor no regime de voo (longe do solo, GPS/baro).
std_sens_pos_voo  = 2.0;            % GPS X,Y,Z     [m]
std_sens_vel_voo  = 0.5;            % GPS Doppler   [m/s]
std_sens_ang_voo  = deg2rad(0.5);   % IMU/AHRS      [rad]
std_sens_taxa_voo = deg2rad(0.1);   % Gyro          [rad/s]

R_ekf_voo = diag([ ...
    std_sens_pos_voo  * ones(1,3), ...
    std_sens_vel_voo  * ones(1,3), ...
    std_sens_ang_voo  * ones(1,3), ...
    std_sens_taxa_voo * ones(1,3), ...
    1e-4 ...   % canal 13 (massa)
].^2);

% ─────────────────────────────────────────────────────────────────────────
% 4. Ruído de Medição R — BAIXA ALTITUDE (RTK + Altímetro de Radar)
% ─────────────────────────────────────────────────────────────────────────
% Mesma estrutura de R_ekf_voo, com 1σ bem menores: perto do solo entram
% RTK e radar, sensores muito mais precisos que o GPS/baro de voo.
std_sens_pos_pouso  = 0.05;           % RTK + radar Z         [m]
std_sens_vel_pouso  = 0.05;           % Radar Doppler         [m/s]
std_sens_ang_pouso  = deg2rad(0.1);   % IMU fundida c/ RTK    [rad]
std_sens_taxa_pouso = deg2rad(0.05);  % Gyro tático           [rad/s]

R_ekf_pouso = diag([ ...
    std_sens_pos_pouso  * ones(1,3), ...
    std_sens_vel_pouso  * ones(1,3), ...
    std_sens_ang_pouso  * ones(1,3), ...
    std_sens_taxa_pouso * ones(1,3), ...
    1e-4 ...
].^2);

% ─────────────────────────────────────────────────────────────────────────
% 5. Altitudes de transição com HISTERESE de 20m
%    O EKF carrega só h_radar_lock pelo workspace (compatibilidade da
%    interface do bloco). A unlock é computada dentro do bloco como
%    h_radar_unlock = h_radar_lock + 20.
% ─────────────────────────────────────────────────────────────────────────
h_radar_lock = 150.0;   % [m]  — entra modo pouso abaixo disso

fprintf('   P0_ekf, Q_ekf carregadas (16 estados).\n');
fprintf('   R_ekf_voo, R_ekf_pouso carregadas (13×13, EKF usa 12×12).\n');
fprintf('   Histerese radar lock: entra a %.0fm, sai a %.0fm.\n', ...
        h_radar_lock, h_radar_lock + 20);
fprintf('   Bias do giroscópio estimado online (estados 14-16).\n');

% ─────────────────────────────────────────────────────────────────────────
% 6. Tracker filter — passa-baixa de 1a ordem no x_hat antes do MPC
% ─────────────────────────────────────────────────────────────────────────
%
%   Constante de tempo τ = 50 ms → corner ≈ 3.2 Hz
%   Justificativa: separação espectral EKF (100 Hz) vs banda do controle (~10 Hz).
%
%   Trade-off de τ:
%     - τ pequeno (5-20 ms)  → pouco atraso, pouco efeito em suavizar
%     - τ médio (30-80 ms)   → balanço, ponto de partida
%     - τ grande (>100 ms)   → suaviza muito, mas come margem de fase
%
tau_tracker = 0.05;   % [s]
dt = 0.001;       % [s] 

fprintf('   tau_tracker carregado: %.3f s (corner ≈ %.1f Hz).\n', ...
        tau_tracker, 1/(2*pi*tau_tracker));