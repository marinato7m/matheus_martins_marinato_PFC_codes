% =========================================================================
% METRICA_QUALIDADE_POUSO.M  —  PLOT + MÉTRICAS PÓS-SIMULAÇÃO
%
%   Script de pós-processamento: rode-o por ÚLTIMO, depois de simular o
%   modelo no Simulink, para avaliar como foi o pouso. Lê a saída 'out' da
%   simulação, imprime as métricas de touchdown no Command Window e gera o
%   painel de gráficos (estado real vs. referência da trajetória).
%
%   PRÉ-REQUISITOS no workspace:
%     out  — SimulationOutput com out.x_sim (estados) e out.u_sim (controle)
%     traj — struct de referência (de gerar_traj_rapido / run_scurve)
%     p    — struct de falcon9_params (limites de empuxo, TVC, etc.)
% =========================================================================

% ── Extração Robusta de Dados (Tempos Independentes) ─────────────────────
% x_sim e u_sim podem ter tempos diferentes (taxas de log distintas), por
% isso cada um guarda seu próprio vetor de tempo.
t_x = out.x_sim.Time; 
t_u = out.u_sim.Time; 

% este bloco normaliza tudo para [amostras × n_estados].
x_raw = out.x_sim.Data;
if ndims(x_raw) == 3
    x = squeeze(x_raw)';
else
    if size(x_raw, 1) == length(t_x)
        x = x_raw;
    else
        x = x_raw';
    end
end

% Mesma normalização para o sinal de controle.
u_raw = out.u_sim.Data;
if ndims(u_raw) == 3
    u_ctrl = squeeze(u_raw)';
else
    if size(u_raw, 1) == length(t_u)
        u_ctrl = u_raw;
    else
        u_ctrl = u_raw';
    end
end

% ── Estado final (touchdown) ─────────────────────────────────────────────
% Última amostra = condição no toque, base de quase todas as métricas.
x_final = x(end, :);
t_final = t_x(end);

fprintf('\n=================================================================\n');
fprintf('  MÉTRICAS DE POUSO\n');
fprintf('=================================================================\n');
fprintf('  Tempo de simulação:     %.2f s\n',   t_final);
fprintf('\n  --- POSIÇÃO ---\n');
fprintf('  x_I (lateral):          %+.3f m\n',  x_final(1));
fprintf('  y_I (lateral):          %+.3f m\n',  x_final(2));
fprintf('  z_I (altitude):         %+.3f m\n',  x_final(3));
% Distância radial ao alvo no plano horizontal (precisão do pouso).
fprintf('  Erro lateral (norm):    %.3f m\n',   sqrt(x_final(1)^2 + x_final(2)^2));
fprintf('\n  --- VELOCIDADE (body) ---\n');
fprintf('  u_B (vx body):          %+.3f m/s\n', x_final(4));
fprintf('  v_B (vy body):          %+.3f m/s\n', x_final(5));
fprintf('  w_B (vz body):          %+.3f m/s\n', x_final(6));   
fprintf('\n  --- ATITUDE ---\n');
fprintf('  phi   (pitch):           %+.2f deg\n', rad2deg(x_final(7)));
fprintf('  theta (yaw):          %+.2f deg\n', rad2deg(x_final(8)));
fprintf('  psi   (roll):            %+.2f deg\n', rad2deg(x_final(9)));
fprintf('\n  --- TAXAS ANGULARES ---\n');
fprintf('  p_B:                    %+.3f deg/s\n', rad2deg(x_final(10)));
fprintf('  q_B:                    %+.3f deg/s\n', rad2deg(x_final(11)));
fprintf('  r_B:                    %+.3f deg/s\n', rad2deg(x_final(12)));
fprintf('\n  --- MASSA ---\n');
fprintf('  Massa final:            %.1f kg\n',   x_final(13));
% Diferença entre massa inicial (p.m0) e final = propelente queimado.
fprintf('  Propelente consumido:   %.1f kg\n',   p.m0 - x_final(13));
fprintf('\n  --- CONTROLE ---\n');
% Picos de atitude ao longo de TODO o voo (não só no toque).
fprintf('  theta_max durante voo:  %.2f deg\n',  rad2deg(max(abs(x(:,8)))));
fprintf('  phi_max durante voo:    %.2f deg\n',  rad2deg(max(abs(x(:,7)))));

% Fração do tempo em que cada eixo de TVC operou perto do batente (>95% do
% curso) — indicador de saturação / falta de autoridade de controle.
sat_dy = sum(abs(u_ctrl(:,2)) > 0.95*p.delta_max) / length(t_u) * 100;
sat_dz = sum(abs(u_ctrl(:,3)) > 0.95*p.delta_max) / length(t_u) * 100;
fprintf('  Saturação TVC delta_y:  %.1f%%\n',   sat_dy);
fprintf('  Saturação TVC delta_x:  %.1f%%\n',   sat_dz);
fprintf('=================================================================\n\n');

% ── Plots ─────────────────────────────────────────────────────────────────
% Painel 2×4: cada subplot compara o estado real (linha cheia) com a
% referência da trajetória (tracejado). Linhas pontilhadas marcam alvos,
% limites de empuxo/TVC e valores de touchdown.
figure('Name','Simulation','NumberTitle','off','Position',[50 50 1400 800]);

% --- Linha 1 ---
% Altitude: real vs. referência.
subplot(2,4,1);
plot(t_x, x(:,3)); hold on;
plot(traj.time, traj.h_d, 'r--');
ylabel('z [m]'); xlabel('t [s]');
legend('real','ref'); title('Altitude');

% Velocidade vertical: real vs. ref, com anotação da velocidade no toque.
subplot(2,4,2);
plot(t_x, x(:,6)); hold on;
plot(traj.time, traj.vz_d, 'r--');
yline(x_final(6), 'k:', sprintf('TD: %.2f m/s', x_final(6)));
ylabel('w_B [m/s]'); xlabel('t [s]');
legend('real','ref'); title('Vertical Velocity');

% Empuxo comandado em kN, com as linhas de T_min e T_max.
subplot(2,4,3);
plot(t_u, u_ctrl(:,1)/1e3); hold on;
yline(p.T_min/1e3, 'r--','T_{min}');
yline(p.T_max/1e3, 'r--','T_{max}');
ylabel('T [kN]'); xlabel('t [s]');
title('Commanded Thrust');

% Velocidades laterais no corpo, com valores de touchdown anotados.
subplot(2,4,4);
plot(t_x, x(:,4), 'b'); hold on;
plot(t_x, x(:,5), 'r');
yline(0, 'k:');
yline(x_final(4), 'b:', sprintf('u_B: %+.2f', x_final(4)));
yline(x_final(5), 'r:', sprintf('v_B: %+.2f', x_final(5)));
ylabel('v [m/s]'); xlabel('t [s]');
legend('u_B (vx)','v_B (vy)');
title('Lateral Velocities');

% --- Linha 2 ---
% Ângulo theta: real vs. referência.
subplot(2,4,5);
plot(t_x, rad2deg(x(:,8))); hold on;
plot(traj.time, rad2deg(traj.theta_d), 'r--');
yline(0,'k:');
ylabel('\theta [deg]'); xlabel('t [s]');
legend('real','ref'); title('Yaw');

% Posições laterais x e y: real vs. referência.
subplot(2,4,6);
plot(t_x, x(:,1), 'b'); hold on;
plot(t_x, x(:,2), 'r');
plot(traj.time, traj.x_d,  'b--');
plot(traj.time, traj.y_d,  'r--');
yline(0,'k:');
ylabel('position [m]'); xlabel('t [s]');
legend('x real','y real','x ref','y ref');
title('Lateral Position');

% Deflexões de TVC nos dois eixos, com os batentes +-delta_max.
subplot(2,4,7);
plot(t_u, rad2deg(u_ctrl(:,2)), 'b'); hold on;
plot(t_u, rad2deg(u_ctrl(:,3)), 'r');
yline( rad2deg(p.delta_max),'k--');
yline(-rad2deg(p.delta_max),'k--');
ylabel('\delta [deg]'); xlabel('t [s]');
legend('\delta_y','\delta_x');
title('TVC Deflection');

% Ângulo phi: real vs. referência.
subplot(2,4,8);
plot(t_x, rad2deg(x(:,7)), 'b'); hold on;
plot(traj.time, rad2deg(traj.phi_d), 'r--');
yline(0,'k:');
ylabel('\phi [deg]'); xlabel('t [s]');
legend('real','ref'); title('Pitch');