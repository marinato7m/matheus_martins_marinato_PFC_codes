% =========================================================================
% TEST_6DOF.M  —  Passo 3: Verificação do Modelo 6-DOF Completo
%
% O QUE O 6-DOF ADICIONA AO 3-DOF
% ─────────────────────────────────
% O 3-DOF validou: cinemática translacional, arrasto, integrador.
% O 6-DOF adiciona:
%   (1) Cinemática rotacional:  euler_dot = L(euler) · omega_B
%   (2) Dinâmica rotacional:    J·omega_dot = M_tvc − omega×(J·omega)
%   (3) Momento do TVC:         M = r_bocal × F_tvc
%   (4) Acoplamento:            v_dot inclui −omega×v_B
%
% CONVENÇÃO DE EIXOS (resumo para os testes)
% ─────────────────────────────────────────────
%   p (x(10)) = taxa angular em torno de X_B  → Mx controla p
%   q (x(11)) = taxa angular em torno de Y_B  → My controla q → afeta theta
%   r (x(12)) = taxa angular em torno de Z_B  → Mz (roll, sem atuador)
%
%   delta_y > 0 → Fy = -T·sin(δy) → Mx < 0 → p < 0 → nariz +Y_B
%   delta_z > 0 → Fx = +T·sin(δz) → My < 0 → q < 0 → nariz +X_B
%
% TESTE E — Conservação de momento angular
%   Valida: termo giroscópico omega×(J·omega) na dinâmica rotacional
%   Referência: Meriam & Kraige (2015), cap. 7
%
% TESTE F — Pulso de torque (alpha = Mx/Jxx)
%   Valida: J, L_TVC e tvc_model produzem momento correto
%   Como: delta_y → Mx = -L·T·sin(δy) → p_dot = Mx/Jxx → p_final = p_dot·t
%   Referência: Tewari (2007), cap. 6 — equação de Euler
%
% TESTE G — Simetria (±delta_y)
%   Valida: convenção de sinais do TVC e do cross(r_bocal, F_tvc)
%   Como: +δy e −δy devem produzir p e x_I espelho exato
%   Referência: simetria física
%
% TESTE H — CG móvel altera autoridade de controle
%   Valida: cg_model comunica L_TVC para dinâmica rotacional
%   Como: mesmo δy, m=m0 vs m=mf → |p| proporcional a L_TVC/Jxx
%   Referência: equação de Euler — sensibilidade a parâmetros
%
% TESTE I — Divergência sem controlador (qualitativo)
%   Valida: modelo reproduz instabilidade open-loop esperada
%   Referência: Tewari (2007), cap. 8
%
% COMO USAR
%   Coloque todos os .m na mesma pasta e execute: test_6dof
% =========================================================================

clear; close all; clc;

p = falcon9_params();

fprintf('\n=========================================================\n');
fprintf('  VERIFICACAO PASSO 3 — MODELO 6-DOF COMPLETO\n');
fprintf('=========================================================\n\n');

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% Estado inicial base: foguete vertical, em repouso
% x = [x_I; y_I; z_I; u_B; v_B; w_B; phi; theta; psi; p_B; q_B; r_B; m]
% Cada teste parte deste estado e ajusta só o que precisa.
x0_base = [0; 0; p.h0;        ... % posição inercial [m]
            0; 0; 0;           ... % velocidade no corpo [m/s]
            0; 0; 0;           ... % ângulos de Euler [rad]
            0; 0; 0;           ... % taxas angulares [rad/s]
            p.m0];

u_off = [0; 0; 0];   % motor desligado, sem deflexão

% =========================================================================
% TESTE E — CONSERVAÇÃO DE MOMENTO ANGULAR
% Valida: omega×(J·omega) na equação de Euler
%
% Sem torques externos (g=0, T=0, sem arrasto), |J·omega| deve ser constante.
% Usamos omega em DOIS eixos (p e q não nulos) para ativar o cruzamento.
% Com um só eixo, omega×(J·omega) = 0 e o teste não valida nada.
% =========================================================================
fprintf('TESTE E — Conservacao de momento angular\n');
fprintf('  Valida: omega x (J*omega) na dinamica rotacional\n');
fprintf('  Referencia: Meriam & Kraige (2015), cap. 7\n');

% Cópia de p sem forças externas, para isolar a dinâmica rotacional livre.
p_E       = p;
p_E.g     = 0;
p_E.CD0   = 0;
p_E.g_I   = [0;0;0];

omega0    = [0.1; 0.05; 0];   % dois eixos ativados — cruzamento vetorial ≠ 0

% Injeta omega inicial em p (x10) e q (x11).
x0_E      = x0_base;
x0_E(10)  = omega0(1);   % p
x0_E(11)  = omega0(2);   % q

[t_E, x_E] = ode45(@(t,x) rocket_6dof(t, x, u_off, p_E), [0, 30], x0_E, opts);

[Jxx_E, Jyy_E, Jzz_E] = inertia_model(p.m0, p_E);

p_E_vec = x_E(:,10);
q_E_vec = x_E(:,11);
r_E_vec = x_E(:,12);

% Módulo do momento angular |J·omega|; deve permanecer constante.
L_mag = sqrt((Jxx_E*p_E_vec).^2 + (Jyy_E*q_E_vec).^2 + (Jzz_E*r_E_vec).^2);
L_ref = L_mag(1);
L_var = max(L_mag) - min(L_mag);
err_E = L_var / L_ref * 100;

fprintf('  |L| inicial:      %.6f kg.m2/s\n', L_ref);
fprintf('  Variacao de |L|:  %.2e kg.m2/s\n', L_var);
fprintf('  Erro relativo:    %.8f%%\n', err_E);
if err_E < 0.01
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% =========================================================================
% TESTE F — PULSO DE TORQUE (alpha = Mx/Jxx)
% Valida: J, L_TVC e tvc_model produzem momento correto
%
% delta_y > 0 → Mx = -L_TVC·T·sin(delta_y) → p_dot = Mx/Jxx (negativo)
% p_final = p_dot · t_pulse
%
% Verificamos o MÓDULO para independer do sinal da convenção.
% O sinal é validado no Teste G (simetria).
% =========================================================================
fprintf('TESTE F — Pulso de torque (alpha = Mx/Jxx)\n');
fprintf('  Valida: J, L_TVC e tvc_model calculam momento corretamente\n');
fprintf('  Referencia: Tewari (2007), cap. 6\n');
fprintf('  Eixo afetado: p_B (x(10)) — delta_y gera Mx que afeta p\n');

p_F     = p;
p_F.g   = 0;
p_F.CD0 = 0;
p_F.g_I = [0;0;0];

delta_y_F = deg2rad(1.0);   % deflexão de 1°
T_F       = p.T_max;
t_pulse   = 1.0;             % s

L_TVC_F      = cg_model(p.m0, p);
[Jxx_F, ~, ~] = inertia_model(p.m0, p);

% Cálculo manual: Mx = -L·T·sin(δy), p_dot = Mx/Jxx
% Previsão analítica de p_final = p_dot·t, comparada com a simulação adiante.
Mx_calc      = -L_TVC_F * T_F * sin(delta_y_F);
p_dot_calc   = Mx_calc / Jxx_F;    % negativo para delta_y > 0
p_final_calc = p_dot_calc * t_pulse;

fprintf('  T = %.0f N,  delta_y = %.4f rad (1 deg)\n', T_F, delta_y_F);
fprintf('  L_TVC = %.2f m,  Jxx = %.3e kg.m2\n', L_TVC_F, Jxx_F);
fprintf('  Mx = -L*T*sin(delta_y) = %.2f N.m  (negativo = correto)\n', Mx_calc);
fprintf('  p_dot = Mx/Jxx = %.6f rad/s2\n', p_dot_calc);
fprintf('  p_final esperado (p_dot*t): %.6f rad/s\n', p_final_calc);

u_F = [T_F; delta_y_F; 0];

[~, x_F] = ode45(@(t,x) rocket_6dof(t, x, u_F, p_F), [0, t_pulse], x0_base, opts);

p_final_sim = x_F(end, 10);   % p_B (índice 10, eixo X_B)

err_F = abs(p_final_sim - p_final_calc) / abs(p_final_calc) * 100;

fprintf('  p_final simulado:           %.6f rad/s\n', p_final_sim);
fprintf('  Sinal: delta_y>0 → p<0 (nariz para +Y_B)?  %s\n', ...
        ternario(p_final_sim < 0, 'SIM (correto)', 'NAO (verificar sinal)'));
fprintf('  Erro modulo:                %.6f%%\n', err_F);
if err_F < 1.0
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% =========================================================================
% TESTE G — SIMETRIA (±delta_y)
% Valida: sinais do TVC e do produto vetorial r_bocal × F_tvc
%
% +δy e −δy devem produzir p_B e x_I espelho exato (soma = 0).
% Teste genuíno agora que sabemos que p_B ≠ 0 (Teste F corrigido).
% =========================================================================
fprintf('TESTE G — Simetria (+-delta_y)\n');
fprintf('  Valida: sinais do TVC e do cross(r_bocal, F_tvc)\n');
fprintf('  Referencia: simetria fisica\n');

p_G     = p;
p_G.g   = 0;
p_G.CD0 = 0;
p_G.g_I = [0;0;0];

delta_G = deg2rad(1.0);
T_G     = p.T_max;
t_G     = 3.0;

% Duas simulações idênticas exceto pelo sinal da deflexão.
u_Gpos = [T_G; +delta_G; 0];
u_Gneg = [T_G; -delta_G; 0];

[t_Gp, x_Gp] = ode45(@(t,x) rocket_6dof(t, x, u_Gpos, p_G), [0, t_G], x0_base, opts);
[t_Gn, x_Gn] = ode45(@(t,x) rocket_6dof(t, x, u_Gneg, p_G), [0, t_G], x0_base, opts);

% Reamostra ambas numa grade comum para comparar ponto a ponto.
t_grid = linspace(0, t_G, 500);
xI_pos = interp1(t_Gp, x_Gp(:,1), t_grid, 'spline');
xI_neg = interp1(t_Gn, x_Gn(:,1), t_grid, 'spline');
p_pos  = interp1(t_Gp, x_Gp(:,10), t_grid, 'spline');   % p_B
p_neg  = interp1(t_Gn, x_Gn(:,10), t_grid, 'spline');

% Simetria: +δy e −δy devem ser espelho (soma ≈ 0)
% A soma das respostas mede o desvio da simetria perfeita.
err_x_abs  = max(abs(xI_pos + xI_neg));
err_p_abs  = max(abs(p_pos  + p_neg));
ref_x      = max(abs(xI_pos));
ref_p      = max(abs(p_pos));
err_x_rel  = err_x_abs / max(ref_x, 1e-10) * 100;
err_p_rel  = err_p_abs / max(ref_p, 1e-10) * 100;

fprintf('  p_B final com +delta: %.6f rad/s\n', p_pos(end));
fprintf('  p_B final com -delta: %.6f rad/s  (esperado: %.6f)\n', p_neg(end), -p_pos(end));
fprintf('  x_I final com +delta: %.4f m\n', xI_pos(end));
fprintf('  x_I final com -delta: %.4f m  (esperado: %.4f)\n', xI_neg(end), -xI_pos(end));
fprintf('  Erro simetria p_B:   %.2e  (%.8f%%)\n', err_p_abs, err_p_rel);
fprintf('  Erro simetria x_I:   %.2e m (%.8f%%)\n', err_x_abs, err_x_rel);
if err_p_rel < 0.01 && err_x_rel < 0.01
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% =========================================================================
% TESTE H — CG MÓVEL ALTERA AUTORIDADE DE CONTROLE
% Valida: cg_model comunica L_TVC corretamente para dinâmica rotacional
%
% Mesmo δy, m=m0 vs m=mf.
% Razão de p_final deve bater com razão (L_TVC_f/Jxx_f) / (L_TVC_0/Jxx_0).
% =========================================================================
fprintf('TESTE H — CG movel altera autoridade de controle\n');
fprintf('  Valida: cg_model afeta L_TVC e portanto p_dot\n');
fprintf('  Eixo verificado: p_B (x(10))\n');

p_H     = p;
p_H.g   = 0;
p_H.CD0 = 0;
p_H.g_I = [0;0;0];

delta_H = deg2rad(1.0);
T_H     = p.T_max;
u_H     = [T_H; delta_H; 0];

% Mesmo comando, dois valores de massa (m0 e mf).
x0_H0        = x0_base;           % m = m0
x0_Hf        = x0_base;
x0_Hf(13)    = p.mf;              % m = mf

[~, x_H0] = ode45(@(t,x) rocket_6dof(t, x, u_H, p_H), [0, 0.5], x0_H0, opts);
[~, x_Hf] = ode45(@(t,x) rocket_6dof(t, x, u_H, p_H), [0, 0.5], x0_Hf, opts);

p_H0 = x_H0(end, 10);   % p_B com m = m0
p_Hf = x_Hf(end, 10);   % p_B com m = mf

% Razão prevista a partir de L_TVC/Jxx em cada massa, vs. razão simulada.
L_TVC_H0 = cg_model(p.m0, p);
L_TVC_Hf = cg_model(p.mf, p);
[Jxx_H0, ~, ~] = inertia_model(p.m0, p);
[Jxx_Hf, ~, ~] = inertia_model(p.mf, p);

razao_prev = (L_TVC_Hf / Jxx_Hf) / (L_TVC_H0 / Jxx_H0);
razao_sim  = p_Hf / p_H0;   % razão das taxas angulares

err_H = abs(razao_sim - razao_prev) / razao_prev * 100;

fprintf('  Em m0=%.0f kg: L_TVC=%.2fm, Jxx=%.3e → p=%.6f rad/s\n', ...
        p.m0, L_TVC_H0, Jxx_H0, p_H0);
fprintf('  Em mf=%.0f kg: L_TVC=%.2fm, Jxx=%.3e → p=%.6f rad/s\n', ...
        p.mf, L_TVC_Hf, Jxx_Hf, p_Hf);
fprintf('  Razao p_f/p_0 simulada:      %.6f\n', razao_sim);
fprintf('  Razao (L_f/J_f)/(L_0/J_0):  %.6f  (previsto)\n', razao_prev);
fprintf('  Erro:                        %.6f%%\n', err_H);
if err_H < 1.0
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% =========================================================================
% TESTE I — DIVERGÊNCIA SEM CONTROLADOR (qualitativo)
% Valida: modelo reproduz instabilidade open-loop esperada
% Um foguete com CG acima do CP deve divergir sem controle.
% =========================================================================
fprintf('TESTE I — Divergencia sem controlador (qualitativo)\n');
fprintf('  Valida: instabilidade open-loop do foguete\n');
fprintf('  Referencia: Tewari (2007), cap. 8\n');

% Pequena perturbação inicial de atitude/taxa, motor sem TVC.
x0_I        = x0_base;
x0_I(8)     = deg2rad(2.0);    % theta0 = 2°
x0_I(11)    = deg2rad(0.5);    % q0 = 0.5°/s

u_I = [p.T_min; 0; 0];   % sem controle (delta = 0)

% Para quando |theta| atinge 45° (divergência confirmada).
opts_I = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @ev_divergencia);

[t_I, x_I, te_I, ~, ~] = ode45(@(t,x) rocket_6dof(t, x, u_I, p), [0, 30], x0_I, opts_I);

theta_I = rad2deg(x_I(:,8));

% Passa se atingiu o limite OU se theta ao menos dobrou (tendência clara).
if ~isempty(te_I)
    fprintf('  Foguete atingiu theta=45 deg em t = %.2f s\n', te_I);
    fprintf('  PASSOU (divergencia confirmada)\n\n');
elseif abs(theta_I(end)) > abs(theta_I(1)) * 2
    fprintf('  theta inicial: %.2f deg\n', theta_I(1));
    fprintf('  theta final:   %.2f deg  (cresceu — divergencia em curso)\n', theta_I(end));
    fprintf('  PASSOU (tendencia de divergencia confirmada)\n\n');
else
    fprintf('  theta NAO divergiu — verificar modelo\n');
    fprintf('  FALHOU\n\n');
end

fprintf('=========================================================\n');
fprintf('  RESUMO\n');
fprintf('  Teste E: |J*omega| conservado (termo girosc. correto)\n');
fprintf('  Teste F: torque correto — p_B bate com Mx/Jxx*t\n');
fprintf('  Teste G: simetria +-delta_y verificada em p_B e x_I\n');
fprintf('  Teste H: cg_model afeta autoridade de controle\n');
fprintf('  Teste I: divergencia open-loop confirmada\n');
fprintf('\n');
fprintf('  PROXIMA ETAPA: implementar controlador.\n');
fprintf('=========================================================\n\n');

% ── PLOTS ─────────────────────────────────────────────────────────────────

% Teste E: taxas angulares (precessão) e |L| constante.
figure('Name','TESTE E — Conservacao de Momento Angular', ...
       'NumberTitle','off','Position',[50 50 800 380]);
subplot(1,2,1);
plot(t_E, x_E(:,10), 'b-', 'LineWidth', 1.5); hold on;
plot(t_E, x_E(:,11), 'g-', 'LineWidth', 1.5);
plot(t_E, x_E(:,12), 'r-', 'LineWidth', 1.5); grid on;
xlabel('Tempo [s]'); ylabel('omega [rad/s]');
legend('p (X_B)','q (Y_B)','r (Z_B)');
title('Taxas angulares (variam por precessao)');
subplot(1,2,2);
plot(t_E, L_mag, 'k-', 'LineWidth', 2); grid on;
yline(L_ref, 'r--', '|L| inicial', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('|J \cdot \omega| [kg\cdotm^2/s]');
title('Teste E — |Momento angular| (deve ser constante)');

% Teste G: respostas a ±δy e suas somas (devem ser ~0 = simetria).
figure('Name','TESTE G — Simetria TVC', ...
       'NumberTitle','off','Position',[50 460 900 380]);
subplot(1,2,1);
plot(t_grid, p_pos,  'b-', 'LineWidth', 2); hold on;
plot(t_grid, p_neg,  'r--','LineWidth', 2);
plot(t_grid, p_pos + p_neg, 'k:', 'LineWidth', 1.5); grid on;
xlabel('Tempo [s]'); ylabel('p_B [rad/s]');
legend('+\delta_y','-\delta_y','soma (deve=0)');
title('Teste G — Simetria em p_B (taxa de roll)');
subplot(1,2,2);
plot(t_grid, xI_pos, 'b-', 'LineWidth', 2); hold on;
plot(t_grid, xI_neg, 'r--','LineWidth', 2);
plot(t_grid, xI_pos + xI_neg, 'k:', 'LineWidth', 1.5); grid on;
xlabel('Tempo [s]'); ylabel('x_I [m]');
legend('+\delta_y','-\delta_y','soma (deve=0)');
title('Teste G — Simetria em x_I');

% Teste I: divergência de pitch e trajetória resultante.
figure('Name','TESTE I — Divergencia Open-Loop', ...
       'NumberTitle','off','Position',[50 870 800 380]);
subplot(1,2,1);
plot(t_I, theta_I, 'b-', 'LineWidth', 2); grid on;
yline(0, 'k--', 'LineWidth', 1);
xlabel('Tempo [s]'); ylabel('\theta [deg]');
title('Teste I — Angulo de Pitch (deve divergir)');
subplot(1,2,2);
plot(t_I, x_I(:,1), 'r-', 'LineWidth', 2); hold on;
plot(t_I, x_I(:,3), 'b-', 'LineWidth', 2); grid on;
xlabel('Tempo [s]'); ylabel('Posicao [m]');
legend('x_I (desvio)','z_I (altitude)');
title('Teste I — Trajetoria divergente');

% ── FUNÇÕES LOCAIS ─────────────────────────────────────────────────────────

function s = ternario(cond, a, b)
% Seleção condicional de valor (operador ternário inexistente em MATLAB).
    if cond; s = a; else; s = b; end
end

function [value, isterminal, direction] = ev_divergencia(~, x)
% Evento de parada: |theta| (x(8)) cruzando 45° de baixo para cima.
    value      = abs(x(8)) - deg2rad(45);
    isterminal = 1;
    direction  = 1;
end