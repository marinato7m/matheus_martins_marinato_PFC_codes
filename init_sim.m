% =========================================================================
% INIT_SIM.M
% Inicializa as variáveis do workspace necessárias para o modelo Simulink.
%
% EXECUTE SEMPRE ANTES DE RODAR A SIMULAÇÃO: primeiro script a rodar, ele
% dá clc e clear all, atenção.
%
% Por que é necessário?
%   Os blocos Simulink referenciam variáveis do workspace MATLAB:
%     - x0_sim : condição inicial do Integrador (13x1)
%     - p      : struct de parâmetros (usada pelos blocos via coder.extrinsic)
%     - p.T_min, p.T_max, etc;
%
% Este é o setup mínimo de malha aberta (open-loop); para a malha com MPC
% usar init_mpc, e para o filtro usar init_ekf depois deste.
% =========================================================================
clc, clear all;
% ── Parâmetros do foguete ──────────────────────────────────────────────────
p = falcon9_params();

% ── Condição inicial do Integrador ────────────────────────────────────────
% Os valores vêm todos do struct p (definidos em falcon9_params), na ordem
% [pos_I; vel_B; euler; taxas_B; massa] esperada pelo Integrador.
x0_sim = [p.x0;              % x_I   [m]
           p.y0;              % y_I   [m]
           p.h0;           % z_I   [m]
           p.vx0;              % u_B   [m/s]
           p.vy0;              % v_B   [m/s]
           p.vz0;          % w_B   [m/s]
           p.phi0;              % phi   [rad]
           p.theta0;   % theta [rad]
           p.psi0;              % psi   [rad]
           p.p0;              % p_B   [rad/s]
           p.q0;   % q_B   [rad/s]
           p.r0;              % r_B   [rad/s]
           p.m0];          % m     [kg]

% ── Variável para o bloco Constant (u_controle) ───────────────────────────
% Open-loop: empuxo mínimo, sem deflexão TVC
% Aqui o comando é fixo (não há controlador na malha): empuxo constante e
% TVC zerado. A linha comentada guarda uma variante com leve deflexão inicial.
% u_openloop = [p.T_min; deg2rad(1); 0];   % [T; delta_y; delta_z]
u_openloop = [p.T_max; 0; 0];   % [T; delta_y; delta_z]

% ── Exibir confirmação ────────────────────────────────────────────────────
fprintf('\n=== Workspace inicializado ===\n');
fprintf('  p         — struct de parametros do foguete\n');
fprintf('  x0_sim    — condicao inicial (13x1)\n');
fprintf('  u_openloop — controle open-loop [%.0fN; 0; 0]\n\n', p.T_min);
fprintf('Estado inicial:\n');
fprintf('  Altitude : %.0f m\n',      x0_sim(3));
fprintf('  Vel. vert: %.1f m/s\n',    x0_sim(6));
fprintf('  Pitch    : %.1f deg\n',    rad2deg(x0_sim(7)));
fprintf('  Yaw    : %.1f deg\n',    rad2deg(x0_sim(8)));
fprintf('  Massa    : %.0f kg\n\n',   x0_sim(13));
fprintf('Rode os init.m restantes e pronto pra simulação no Simulink\n');
