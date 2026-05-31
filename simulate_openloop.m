% =========================================================================
% SIMULATE_OPENLOOP.M
% Simula o foguete em malha aberta (sem controlador) com o modelo 6-DOF.
%
% OBJETIVO
%   Ter uma base visual de como o foguete se comporta SEM controle.
%   Usado para:
%     - Ter referência de comparação quando o controlador for implementado
%     - Debuggar problemas no modelo antes de adicionar o controlador
%
% ENTRADAS (configuráveis na Seção 1)
%   T_const  — empuxo constante durante a simulação [N]
%   delta_y0 — deflexão TVC constante em Y [rad]  (0 = sem TVC)
%   delta_z0 — deflexão TVC constante em Z [rad]  (0 = sem TVC)
%   theta0   — inclinação inicial [deg]  (0 = vertical perfeito)
%
% SAÍDA
%   t, x, u_hist — ficam no workspace para inspeção
%   Figuras 1-4  — geradas por plot_6dof()
%
% COMO USAR
%   Coloque todos os .m na mesma pasta e execute: simulate_openloop
% =========================================================================

clear; close all; clc;

p = falcon9_params();

% =========================================================================
%  SEÇÃO 1 — CONFIGURAÇÃO  (edite aqui)
% =========================================================================

% Empuxo constante durante a simulação; escolha, 0, T_min ou T_max
T_const = p.T_min;      % [N]

% Deflexão TVC constante (0 = sem controle de atitude)
delta_y0 = deg2rad(0);  % [rad]
delta_z0 = deg2rad(0);  % [rad]

% Condições iniciais de atitude
theta0_deg = 2.0;       % inclinação inicial de pitch [graus]
q0_dps     = 0.5;       % taxa inicial de pitch [graus/s]

% Tempo máximo de simulação
% O ode45 para antes se: (a) z_I < 0 (bateu no solo) ou (b) |theta| > 45°
t_max = 30.0;           % [s]

% =========================================================================
%  SEÇÃO 2 — ESTADO INICIAL
% =========================================================================

x0 = [0;                   ... % x_I [m]
       0;                   ... % y_I [m]
       p.h0;               ... % z_I [m]
       0;                   ... % u_B [m/s]
       0;                   ... % v_B [m/s]
       p.vz0;              ... % w_B [m/s]  — velocidade vertical inicial
       0;                   ... % phi [rad]
       deg2rad(theta0_deg); ... % theta [rad]
       0;                   ... % psi [rad]
       0;                   ... % p_B [rad/s]
       deg2rad(q0_dps);    ... % q_B [rad/s]
       0;                   ... % r_B [rad/s]
       p.m0];                   % m [kg]

% Entrada constante ao longo da simulação; verificar se é o valor dentro do
% bloco Const no simulink
u_const = [T_const; delta_y0; delta_z0];

% =========================================================================
%  SEÇÃO 3 — INTEGRAÇÃO
% =========================================================================

fprintf('\n=========================================================\n');
fprintf('  SIMULAÇÃO OPEN-LOOP — 6-DOF\n');
fprintf('  T = %.0f kN (%.0f%% throttle)\n', T_const/1e3, T_const/p.T_max*100);
fprintf('  theta0 = %.1f deg,  q0 = %.1f deg/s\n', theta0_deg, q0_dps);
fprintf('  Sem controlador — divergência esperada\n');
fprintf('=========================================================\n\n');

% Tolerâncias apertadas e eventos de parada (solo / divergência de atitude).
opts = odeset('RelTol',   1e-8, ...
              'AbsTol',   1e-8, ...
              'Events',   @(t,x) eventos_parada(t, x, p));

% Integra o modelo 6-DOF com entrada constante; mede o tempo de cálculo.
tic;
[t, x, te, xe, ie] = ode45(@(t,x) rocket_6dof(t, x, u_const, p), ...
                             [0, t_max], x0, opts);
tempo_calc = toc;

% =========================================================================
%  SEÇÃO 4 — RELATÓRIO DO CONSOLE
% =========================================================================

fprintf('Simulação concluída em %.2f s (tempo de cálculo)\n\n', tempo_calc);

% Informa qual evento parou a simulação (ou se rodou até t_max).
if ~isempty(te)
    motivos = {'Solo (z_I < 0)', 'Divergência (|theta| > 45°)'};
    fprintf('Simulação parou em t = %.2f s\n', te(1));
    fprintf('Motivo: %s\n\n', motivos{ie(1)});
else
    fprintf('Simulação completou t_max = %.1f s sem evento de parada\n\n', t_max);
end

% Estado final
fprintf('─── Estado Final ────────────────────────────────────────\n');
fprintf('  Altitude:          %.2f m\n',   x(end,3));
fprintf('  Velocidade vert.:  %.2f m/s\n', x(end,6));
fprintf('  theta:       %.2f deg\n', rad2deg(x(end,8)));
fprintf('  q:      %.2f deg/s\n', rad2deg(x(end,11)));
fprintf('  Massa final:       %.1f kg\n',  x(end,13));
fprintf('  Propelente usado:  %.1f kg\n',  p.m0 - x(end,13));
fprintf('  Desvio lateral x:  %.2f m\n',  x(end,1));
fprintf('─────────────────────────────────────────────────────────\n\n');

% =========================================================================
%  SEÇÃO 5 — HISTÓRICO DE ENTRADAS (para plot_6dof)
% =========================================================================

% u_hist é constante nesta simulação open-loop
% Replica u_const em todas as N amostras de tempo (formato esperado por plot_6dof).
u_hist = repmat(u_const', length(t), 1);   % N×3

% =========================================================================
%  SEÇÃO 6 — PLOTS
% =========================================================================

nome_sim = sprintf('Open-loop (T=%.0fkN, theta0=%.1fdeg)', ...
                   T_const/1e3, theta0_deg);

figs = plot_6dof(t, x, u_hist, p, nome_sim);

% ── Plot adicional: theta(t) em destaque (o mais importante no open-loop)
% Painel dedicado à divergência de pitch e ao desvio horizontal que ela causa.
figure('Name', 'Divergência de Pitch — Open-Loop', ...
       'NumberTitle', 'off', 'Position', [700 50 600 450]);

% Pitch com as linhas de vertical perfeito e do limite de ±45°.
subplot(2,1,1);
plot(t, rad2deg(x(:,8)), 'Color', [0.49 0.18 0.56], 'LineWidth', 2.5);
hold on;
yline(0,  'k--', 'Vertical perfeito', 'LineWidth', 1);
yline(45, 'r--', '45° — limite', 'LineWidth', 1.5);
yline(-45,'r--', '', 'LineWidth', 1.5);
xlabel('Tempo [s]'); ylabel('\theta [deg]');
title('Ângulo de Pitch — Divergência Open-Loop');
grid on;

% Desvio horizontal: consequência direta da inclinação não corrigida.
subplot(2,1,2);
plot(t, x(:,1), 'Color', [0.85 0.33 0.10], 'LineWidth', 2); hold on;
plot(t, x(:,2), 'Color', [0.00 0.45 0.70], 'LineWidth', 2);
legend('x_I (Norte)', 'y_I (Leste)', 'Location', 'best');
xlabel('Tempo [s]'); ylabel('Posição [m]');
title('Desvio Horizontal — Causado pela Inclinação');
grid on;

sgtitle('DIVERGÊNCIA OPEN-LOOP — Sem Controlador', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');

% =========================================================================
%  SEÇÃO 7 — EXPORTAR PARA WORKSPACE
% =========================================================================

% Variáveis disponíveis no workspace após a execução:
%   t      — vetor de tempo [s]
%   x      — matriz de estados (N×13)
%   u_hist — entradas (N×3)
%   p      — parâmetros do foguete
%   te, xe — tempo e estado do evento de parada (se houver)

fprintf('Variáveis no workspace: t, x, u_hist, p\n');
fprintf('Para plotar novamente: plot_6dof(t, x, u_hist, p, nome_sim)\n\n');

% ── FUNÇÕES LOCAIS ─────────────────────────────────────────────────────────

function [value, isterminal, direction] = eventos_parada(~, x, p)
% Dois eventos de parada:
%   Evento 1: z_I < 0  (foguete bateu no solo)
%   Evento 2: |theta| > 45°  (divergência severa — simulação inútil a partir daqui)
% Cada componente de 'value' cruza zero quando o evento ocorre.

    z_I   = x(3);
    theta = x(8);

    value      = [z_I;                        ... % evento 1: z cruzando zero
                  abs(theta) - deg2rad(45)];       % evento 2: |theta| cruzando 45°
    isterminal = [1; 1];    % ambos param a simulação
    direction  = [-1; 1];   % z decrescendo; |theta| crescendo
end