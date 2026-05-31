function figs = plot_6dof(t, x, u_hist, p, nome, ref)
% Esta função é chamada por simulate_openloop.m e, se  quiser, pelo
% controlador. O formato de saída é idêntico — facilita comparação
% direta entre "sem controle" e "com controle" no mesmo estilo de plot.
%
% ENTRADAS
%   t      — vetor de tempo [s]  (N×1)
%   x      — matriz de estados (N×13)
%              x(:,1:3)  = [x_I, y_I, z_I]    posição inercial [m]
%              x(:,4:6)  = [u_B, v_B, w_B]    velocidade no corpo [m/s]
%              x(:,7:9)  = [phi, theta, psi]  ângulos de Euler [rad]
%              x(:,10:12)= [p_B, q_B, r_B]   taxas angulares [rad/s]
%              x(:,13)   = m                  massa [kg]
%   u_hist — entradas de controle ao longo do tempo (N×3)
%              u_hist(:,1) = T [N]
%              u_hist(:,2) = delta_y [rad]
%              u_hist(:,3) = delta_z [rad]
%   p      — struct de parâmetros (falcon9_params)
%   nome   — string com nome da simulação (aparece nos títulos)
%   ref    — (opcional) struct com trajetória de referência S-curve:
%              ref.time, ref.h_d, ref.vz_d, ref.az_d
%
% SAÍDA
%   figs — vetor com handles das figuras geradas

% ── Extrair estados ────────────────────────────────────────────────────────
x_I   = x(:,1);   y_I   = x(:,2);   z_I   = x(:,3);
u_B   = x(:,4);   v_B   = x(:,5);   w_B   = x(:,6);
phi   = rad2deg(x(:,7));
theta = rad2deg(x(:,8));
psi   = rad2deg(x(:,9));
p_B   = rad2deg(x(:,10));
q_B   = rad2deg(x(:,11));
r_B   = rad2deg(x(:,12));
m     = x(:,13);

% ── Grandezas derivadas ────────────────────────────────────────────────────
% Velocidade vertical no inercial
vz_I = zeros(size(t));   % velocidade vertical no inercial
for i = 1:length(t)
    euler_i = x(i,7:9)';
    vel_Bi  = x(i,4:6)';
    R = rotation_matrix(euler_i(1), euler_i(2), euler_i(3));
    v_I_i = R * vel_Bi;
    vz_I(i) = v_I_i(3);
end

v_mag      = sqrt(u_B.^2 + v_B.^2 + w_B.^2);   % velocidade [m/s]
prop_usado = p.m0 - m;                            % propelente consumido [kg]
% L_TVC e inércias reconstruídos a cada instante a partir da massa, para
% mostrar como variam ao longo do burn.
L_TVC_t    = arrayfun(@(mi) cg_model(mi, p), m);  % L_TVC ao longo do tempo
[Jxx_t, ~, Jzz_t] = arrayfun(@(mi) inertia_model(mi, p), m);

% ── Referência disponível? ─────────────────────────────────────────────────
% Só sobrepõe a S-curve se o argumento ref foi passado e não está vazio.
tem_ref = nargin >= 6 && ~isempty(ref);

% ── Cores e estilo ─────────────────────────────────────────────────────────
% Paleta fixa por grandeza, para que as 4 figuras fiquem visualmente
% consistentes entre si e entre simulações diferentes.
C.az  = [0.00 0.45 0.70];   % azul escuro — altitude/z
C.vz  = [0.85 0.33 0.10];   % laranja — velocidade
C.phi = [0.47 0.67 0.19];   % verde — roll
C.th  = [0.49 0.18 0.56];   % roxo — pitch
C.ps  = [0.30 0.75 0.93];   % ciano — yaw
C.ref = [0.90 0.00 0.00];   % vermelho — referência
C.ctrl= [0.50 0.50 0.50];   % cinza — entradas

% Helper que prefixa todos os títulos com o nome da simulação.
titulo = @(s) sprintf('[%s]  %s', nome, s);
figs   = gobjects(4,1);

% =========================================================================
% FIGURA 1 — TRANSLAÇÃO: altitude, velocidade vertical, trajetória
% =========================================================================
figs(1) = figure('Name', sprintf('%s — Traducao', nome), ...
                 'NumberTitle', 'off', ...
                 'Position', [50 50 1000 700]);

% Altitude (com solo marcado e, se houver, referência S-curve).
subplot(3,2,1);
plot(t, z_I, '-', 'Color', C.az, 'LineWidth', 2); hold on;
if tem_ref
    plot(ref.time, ref.h_d, '--', 'Color', C.ref, 'LineWidth', 1.5);
    legend('Simulado', 'S-curve ref.', 'Location', 'best');
end
yline(0, 'k--', 'Solo', 'LineWidth', 1);
xlabel('Tempo [s]'); ylabel('z_I [m]');
title(titulo('Altitude')); grid on;

% Velocidade vertical inercial vs. referência.
subplot(3,2,2);
plot(t, vz_I, '-', 'Color', C.vz, 'LineWidth', 2); hold on;
if tem_ref
    plot(ref.time, ref.vz_d, '--', 'Color', C.ref, 'LineWidth', 1.5);
    legend('Simulado', 'S-curve ref.', 'Location', 'best');
end
yline(0, 'k--', 'v=0', 'LineWidth', 1);
xlabel('Tempo [s]'); ylabel('vz [m/s]');
title(titulo('Velocidade Vertical (inercial)')); grid on;

% Desvio horizontal nos dois eixos inerciais.
subplot(3,2,3);
plot(t, x_I, '-', 'Color', C.az, 'LineWidth', 2); hold on;
plot(t, y_I, '--', 'Color', C.vz, 'LineWidth', 2);
legend('x_I (Norte)', 'y_I (Leste)', 'Location', 'best');
xlabel('Tempo [s]'); ylabel('Posição [m]');
title(titulo('Desvio Horizontal')); grid on;

% Módulo da velocidade (norma das 3 componentes do corpo).
subplot(3,2,4);
plot(t, v_mag, '-', 'Color', C.az, 'LineWidth', 2);
xlabel('Tempo [s]'); ylabel('|v| [m/s]');
title(titulo('Speed (norma da velocidade)')); grid on;

% Massa, com a massa final de referência.
subplot(3,2,5);
plot(t, m, '-', 'Color', C.az, 'LineWidth', 2); hold on;
yline(p.mf, 'r--', sprintf('m_f = %.0f kg', p.mf), 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('Massa [kg]');
title(titulo('Massa')); grid on;

% Propelente consumido (m0 - m).
subplot(3,2,6);
plot(t, prop_usado, '-', 'Color', C.vz, 'LineWidth', 2);
xlabel('Tempo [s]'); ylabel('Propelente [kg]');
title(titulo('Propelente Consumido')); grid on;

sgtitle(sprintf('TRANSLAÇÃO — %s', nome), 'FontSize', 12, 'FontWeight', 'bold');

% =========================================================================
% FIGURA 2 — ROTAÇÃO: ângulos de Euler e taxas angulares
% =========================================================================
% Linha de cima: ângulos (roll/pitch/yaw); linha de baixo: suas taxas.
figs(2) = figure('Name', sprintf('%s — Rotacao', nome), ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 1000 600]);

subplot(2,3,1);
plot(t, phi, '-', 'Color', C.phi, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Tempo [s]'); ylabel('\phi [deg]');
title(titulo('Pitch')); grid on;

subplot(2,3,2);
plot(t, theta, '-', 'Color', C.th, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Tempo [s]'); ylabel('\theta [deg]');
title(titulo('Yaw')); grid on;

subplot(2,3,3);
plot(t, psi, '-', 'Color', C.ps, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Tempo [s]'); ylabel('\psi [deg]');
title(titulo('Roll')); grid on;

subplot(2,3,4);
plot(t, p_B, '-', 'Color', C.phi, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Tempo [s]'); ylabel('p [deg/s]');
title(titulo('Taxa de Pitch')); grid on;

subplot(2,3,5);
plot(t, q_B, '-', 'Color', C.th, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Tempo [s]'); ylabel('q [deg/s]');
title(titulo('Taxa de Yaw')); grid on;

subplot(2,3,6);
plot(t, r_B, '-', 'Color', C.ps, 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Tempo [s]'); ylabel('r [deg/s]');
title(titulo('Taxa de Roll')); grid on;

sgtitle(sprintf('ROTAÇÃO — %s', nome), 'FontSize', 12, 'FontWeight', 'bold');

% =========================================================================
% FIGURA 3 — ENTRADAS DE CONTROLE
% =========================================================================
% Só gera se houver histórico de entradas; mostra empuxo e as duas
% deflexões de TVC com seus limites.
if ~isempty(u_hist)
    figs(3) = figure('Name', sprintf('%s — Entradas', nome), ...
                     'NumberTitle', 'off', ...
                     'Position', [150 150 1000 450]);

    % Empuxo em kN, com T_min e T_max.
    subplot(1,3,1);
    plot(t, u_hist(:,1)/1e3, '-', 'Color', C.az, 'LineWidth', 2); hold on;
    yline(p.T_max/1e3, 'r--', 'T_{max}', 'LineWidth', 1.2);
    yline(p.T_min/1e3, 'b--', 'T_{min}', 'LineWidth', 1.2);
    xlabel('Tempo [s]'); ylabel('T [kN]');
    title(titulo('Empuxo')); grid on;

    % Deflexão TVC em Y, com batentes ±delta_max.
    subplot(1,3,2);
    plot(t, rad2deg(u_hist(:,2)), '-', 'Color', C.th, 'LineWidth', 2); hold on;
    yline( rad2deg(p.delta_max), 'r--', '\delta_{max}', 'LineWidth', 1.2);
    yline(-rad2deg(p.delta_max), 'r--', '', 'LineWidth', 1.2);
    xlabel('Tempo [s]'); ylabel('\delta_y [deg]');
    title(titulo('Deflexão TVC — Y')); grid on;

    % Deflexão TVC no outro eixo, com os mesmos batentes.
    subplot(1,3,3);
    plot(t, rad2deg(u_hist(:,3)), '-', 'Color', C.ps, 'LineWidth', 2); hold on;
    yline( rad2deg(p.delta_max), 'r--', '\delta_{max}', 'LineWidth', 1.2);
    yline(-rad2deg(p.delta_max), 'r--', '', 'LineWidth', 1.2);
    xlabel('Tempo [s]'); ylabel('\delta_x [deg]');
    title(titulo('Deflexão TVC — X')); grid on;

    sgtitle(sprintf('ENTRADAS DE CONTROLE — %s', nome), ...
            'FontSize', 12, 'FontWeight', 'bold');
end

% =========================================================================
% FIGURA 4 — PARÂMETROS VARIÁVEIS (CG, inércia)
% =========================================================================
% Mostra como o braço de momento do TVC e as inércias mudam com a massa —
% explica a perda de autoridade de controle ao final do burn.
figs(4) = figure('Name', sprintf('%s — Parametros Variaveis', nome), ...
                 'NumberTitle', 'off', ...
                 'Position', [200 200 900 400]);

% Braço de momento do TVC (= distância CG-bocal), com os extremos r_CG.
subplot(1,3,1);
plot(t, L_TVC_t, '-', 'Color', C.az, 'LineWidth', 2); hold on;
yline(p.r_CG_0, 'r--', sprintf('r_{CG,0}=%.1fm', p.r_CG_0), 'LineWidth', 1.2);
yline(p.r_CG_f, 'b--', sprintf('r_{CG,f}=%.1fm', p.r_CG_f), 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('L_{TVC} [m]');
title(titulo('Braço de Momento TVC')); grid on;

% Inércia de pitch (Jxx), normalizada para 10^6.
subplot(1,3,2);
plot(t, Jxx_t/1e6, '-', 'Color', C.th, 'LineWidth', 2);
xlabel('Tempo [s]'); ylabel('J_{xx} [10^6 kg·m^2]');
title(titulo('Inércia de Pitch (Jxx)')); grid on;

% Inércia de roll (Jzz), normalizada para 10^4.
subplot(1,3,3);
plot(t, Jzz_t/1e4, '-', 'Color', C.phi, 'LineWidth', 2);
xlabel('Tempo [s]'); ylabel('J_{zz} [10^4 kg·m^2]');
title(titulo('Inércia de Roll (Jzz)')); grid on;

sgtitle(sprintf('PARÂMETROS VARIÁVEIS COM MASSA — %s', nome), ...
        'FontSize', 12, 'FontWeight', 'bold');

fprintf('[plot_6dof] %s — 4 figuras geradas.\n', nome);
end