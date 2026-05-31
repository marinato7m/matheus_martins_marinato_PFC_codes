% =========================================================================
% INIT_MPC.M  —  Inicialização do MPC para o modelo Simulink
%
% Execute DEPOIS de init_sim.m
%
% O que este script faz:
%   1. Carrega parâmetros do Falcon 9
%   2. Define Q_mpc e R_mpc
%   3. Define x0_sim (redundância mesmo)
%   4. Cria opts_qp (opções do quadprog — NECESSÁRIO pro MPC)
% =========================================================================

p = falcon9_params();

% ── Condição inicial da simulação ───────────────────────────────────────
% redundancia
x0_sim = [p.x0; p.y0; p.h0;
          p.vx0; p.vy0; p.vz0;
          p.phi0; p.theta0; p.psi0;
          p.p0; p.q0; p.r0;
          p.m0];

% ── Pesos Q e R ──────────────
%
%   Regra de Bryson: Q_ii = 1/desvio_max²
%
% ── Multiplicadores finais, vindos da otimização 2-estágios (LHS+bayesopt).
%    Cada knob escala um grupo de termos de Q/R sobre a referência de Bryson.
k_pos    = 41.8856;
k_vel    = 0.4940;
k_h      = 11.2475;
k_vz     = 1.0737;
k_ang    = 0.4550;
k_rate   = 0.2210;
k_tvc    = 0.2170;

% Matriz de estado Q (10×10): peso = knob / (desvio_max)^2 por estado.
% Posição e velocidade laterais repetem o mesmo knob (simetria x/y).
  Q_mpc = diag([...
      k_pos  / (5.0)^2,           ... %  1. x_I
      k_pos  / (5.0)^2,           ... %  2. y_I
      k_h    / (2.0)^2,           ... %  3. z_I
      k_vel  / (3.0)^2,           ... %  4. u_B 3
      k_vel  / (3.0)^2,           ... %  5. v_B 3
      k_vz   / (1.0)^2,           ... %  6. w_B
      k_ang  / (deg2rad(5))^2,    ... %  7. phi
      k_ang  / (deg2rad(5))^2,    ... %  8. theta
      k_rate / (deg2rad(5))^2,    ... %  9. p_B
      k_rate / (deg2rad(5))^2     ... % 10. q_B
  ]);

% Matriz de controle R (3×3): empuxo com peso base fixo; as duas deflexões
% de TVC compartilham k_tvc, normalizadas pelo curso máximo do gimbal.
  R_mpc = diag([...
      1/(p.T_max - p.T_min)^2,    ... % T
      k_tvc / (p.delta_max)^2,     ... % delta_y
      k_tvc / (p.delta_max)^2      ... % delta_z
  ]);

  % ── Rate weighting (peso S — penaliza variação do comando) ──────────────
%
%   "input rate penalty", padrão industrial em MPC.
%   Penaliza δu = u(k) - u(k-1) para evitar chattering induzido por ruído
%   de estimação
k_rate_T     = 300;     % multiplicador rate do empuxo
k_rate_delta = 50;    % multiplicador rate do TVC 

% Normalização análoga à de R, mas aplicada à taxa de variação do comando.
S_rate = diag([ ...
    k_rate_T     / (p.T_max - p.T_min)^2, ...    % d(T)/dt
    k_rate_delta / (p.delta_max)^2,        ...    % d(δy)/dt
    k_rate_delta / (p.delta_max)^2          ...    % d(δz)/dt
]);


% ── Opções do quadprog (necessário pro MPC) ─────────────────────
%
%   O controlador_mpc.m carrega isto via evalin('base', 'opts_qp').
%   Sem esta variável no workspace, o MPC vai dar erro.
%
%   MaxIterations limita o custo no pior caso (o QP roda a cada passo);
%   interior-point-convex é adequado ao QP convexo do MPC.
%
opts_qp = optimoptions('quadprog', ...
    'Display',   'off', ...
    'Algorithm', 'interior-point-convex', ...
    'MaxIterations', 50, ...           % limita tempo no worst-case
    'OptimalityTolerance', 1e-6);

fprintf('\n=========================================================\n');
fprintf('  INIT_MPC — Workspace pronto para simulação\n');
fprintf('=========================================================\n');
fprintf('  x0_sim:  [%.1f, %.1f, %.1f] m  |  vz=%.1f m/s\n', ...
        x0_sim(1), x0_sim(2), x0_sim(3), x0_sim(6));
fprintf('  Atitude:  phi=%.1f°  theta=%.1f°\n', ...
        rad2deg(x0_sim(7)), rad2deg(x0_sim(8)));
fprintf('  Pesos:   k_vel=%.5f  k_pos=%.5f  k_h=%.5f  k_tvc=%.5f k_vz=%.5f k_ang=%.5f k_rate=%.5f\n', ...
        k_vel, k_pos, k_h, k_tvc, k_vz, k_ang, k_rate);
fprintf('  opts_qp: interior-point-convex, MaxIter=50\n');
fprintf('=========================================================\n\n');
fprintf('Rode o init_ekf.m e pronto pra simulação no Simulink\n');
