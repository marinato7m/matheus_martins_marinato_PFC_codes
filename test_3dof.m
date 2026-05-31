% =========================================================================
% TEST_3DOF.M  —  Passo 2: Verificação do Modelo 3-DOF Translacional
%
% O QUE O 3-DOF ADICIONA AO 1-DOF
% ─────────────────────────────────
% O 1-DOF validou: gravidade, empuxo, consumo de massa, integrador.
% O 3-DOF adiciona exatamente duas coisas novas:
%   (1) Movimento em 3 eixos — velocidade lateral (u_B, v_B) se integra
%       em posição inercial (x_I, y_I) via ṗ_I = R_IB · v_B
%   (2) Arrasto 3D — |v_rel|² agora usa as 3 componentes de velocidade
%
% DISTINÇÃO IMPORTANTE — 3-DOF vs 6-DOF
% ───────────────────────────────────────
% No 3-DOF: atitude é PARÂMETRO FIXO — R_IB constante, não evolui.
% No 6-DOF: atitude é ESTADO DINÂMICO — R_IB muda porque os ângulos de
%           Euler evoluem pela equação cinemática de rotação.
%
% Por isso todos os testes aqui usam theta=0 (R_IB = eye(3), foguete
% vertical). Assim testamos apenas o que é novo no 3-DOF sem ambiguidade
% com o 6-DOF. Testes com atitude não-nula ficam para o test_6dof.m.
%
% TESTE B — Projétil com velocidade lateral (theta=0)
%   Valida: velocidade lateral u_B integra corretamente em x_I
%   Com theta=0: R_IB = eye(3), então ṗ_I = v_B diretamente.
%   Solução analítica: x(t) = u0·t,  z(t) = z0 + vz0·t − ½g·t²
%   Referência: Meriam & Kraige (2015), cap. 2.4
%   Novidade testada: existência do eixo x_I e integração de u_B → x_I
%
% TESTE C — Conservação de energia em translação 3D
%   Valida: modelo não cria nem dissipa energia com movimento lateral
%   E = ½(u² + v² + w²) + g·z deve ser constante (sem empuxo, sem arrasto)
%   Referência: Tewari (2007)
%   Novidade testada: |v|² agora tem 3 componentes
%
% TESTE D — Velocidade terminal com arrasto
%   Valida: aero_model.m computa D = ½ρv²·Sref·CD0 corretamente
%   Na velocidade terminal: D = m·g → az = 0 → vz = constante
%   v_t = √(2mg / ρ·Sref·CD0)
%   Referência: Anderson (2015), cap. 5; Sutton & Biblarz (2017), cap. 2
%   Novidade testada: o modelo de arrasto (não existia no 1-DOF)
%
% COMO USAR
%   Coloque todos os .m na mesma pasta e execute: test_3dof
% =========================================================================

clear; close all; clc;

p = falcon9_params();

fprintf('\n=========================================================\n');
fprintf('  VERIFICACAO PASSO 2 — MODELO 3-DOF TRANSLACIONAL\n');
fprintf('  Todos os testes com theta=0 (R_IB=eye(3), sem ambiguidade\n');
fprintf('  com 6-DOF). Atitude nao-nula fica para test_6dof.m.\n');
fprintf('=========================================================\n\n');

% ODE 3-DOF com atitude FIXA theta=0
% Atitude entra como parâmetro (não estado); com euler nulo, R_IB = I.
euler_fixo = [0; 0; 0];   % foguete vertical, R_IB = eye(3)
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% =========================================================================
% TESTE B — PROJÉTIL COM VELOCIDADE LATERAL
% Valida: u_B (velocidade lateral no corpo) integra corretamente em x_I
%
% Setup:
%   theta=0 → R_IB = eye(3) → ṗ_I = v_B diretamente (sem rotação)
%   Velocidade inicial: u0 = 50 m/s lateral, w0 = −50 m/s vertical
%
% Solução analítica (Meriam & Kraige, 2015, cap. 2.4):
%   x(t) = u0·t
%   z(t) = z0 + w0·t − ½g·t²
%
% Com theta=0, u_B = vx_I exatamente — não há decomposição por R_IB.
% O que testamos: que o integrador cria e popula corretamente o eixo x_I,
% e que a cinemática ṗ_I = R_IB·v_B funciona em 3D.
% =========================================================================
fprintf('TESTE B — Projetil com velocidade lateral (theta=0)\n');
fprintf('  Valida: u_B integra em x_I; cinemática 3D sem rotação\n');
fprintf('  Referencia: Meriam & Kraige (2015), cap. 2.4\n');
fprintf('  Nota: R_IB=eye(3) aqui — teste de 3D puro, sem ambiguidade 6-DOF\n');

u0  =  50.0;    % velocidade lateral inicial [m/s]  → gera x_I
w0  = -50.0;    % velocidade vertical inicial [m/s] (negativo = descendo)
z0  = p.h0;     % altitude inicial [m]

% Tempo de impacto (quadrática em z):
t_hit = (w0 + sqrt(w0^2 + 2*p.g*z0)) / p.g;

% Soluções analíticas de projétil para comparação.
t_vec  = linspace(0, t_hit, 1000);
x_anal = u0 * t_vec;                                    % x analítico
z_anal = z0 + w0*t_vec - 0.5*p.g*t_vec.^2;             % z analítico

% Sem arrasto para este teste (CD0=0)
% Cópia de p com arrasto desligado, para isolar a cinemática pura.
p_nodrag     = p;
p_nodrag.CD0 = 0;

% Estado inicial: [x;y;z; u_B;v_B;w_B; m]
% Com theta=0: u_B=50 m/s lateral e w_B=−50 m/s vertical
x0_B = [0; 0; z0;     u0; 0; w0;     p.m0];

[t_B, x_B] = ode45(@(t,x) ode_3dof(t, x, euler_fixo, 0, p_nodrag), ...
                    [0, t_hit], x0_B, opts);

% Reamostra numérico na grade analítica e mede erros relativos.
xI_num = interp1(t_B, x_B(:,1), t_vec, 'spline');
zI_num = interp1(t_B, x_B(:,3), t_vec, 'spline');

err_x = max(abs(xI_num - x_anal)) / max(abs(x_anal)) * 100;
err_z = max(abs(zI_num - z_anal)) / abs(z0) * 100;

fprintf('  t_hit = %.4f s\n', t_hit);
fprintf('  Alcance horizontal esperado: %.2f m\n', u0 * t_hit);
fprintf('  Alcance horizontal simulado: %.2f m\n', x_B(end,1));
fprintf('  Erro max. x_I: %.2e m  (%.8f%%)\n', max(abs(xI_num-x_anal)), err_x);
fprintf('  Erro max. z_I: %.2e m  (%.8f%%)\n', max(abs(zI_num-z_anal)), err_z);
if err_x < 0.01 && err_z < 0.01
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% =========================================================================
% TESTE C — CONSERVAÇÃO DE ENERGIA EM TRANSLAÇÃO 3D
% Valida: energia mecânica conservada com movimento lateral
%
% Usa os dados do Teste B (sem arrasto, sem empuxo).
% E = ½(u² + v² + w²) + g·z deve ser constante.
% Novidade em relação ao Teste 4 do 1-DOF: |v|² tem componente lateral.
% =========================================================================
fprintf('TESTE C — Conservacao de energia em 3D\n');
fprintf('  Valida: energia conservada com movimento lateral (u_B nao zero)\n');
fprintf('  Referencia: Tewari (2007)\n');

% Reaproveita a trajetória do Teste B; energia inclui as 3 componentes.
u_B   = x_B(:,4);
v_B_c = x_B(:,5);
w_B   = x_B(:,6);
z_I   = x_B(:,3);

e_num = 0.5*(u_B.^2 + v_B_c.^2 + w_B.^2) + p.g*z_I;
e_var = max(e_num) - min(e_num);
e_ref = abs(0.5*(u0^2 + w0^2) + p.g*z0);
err_C = e_var / e_ref * 100;

fprintf('  |v|^2 inicial = u0^2 + w0^2 = %.0f + %.0f = %.0f m2/s2\n', ...
        u0^2, w0^2, u0^2+w0^2);
fprintf('  Variacao de energia: %.2e J/kg\n', e_var);
fprintf('  Erro relativo:       %.10f%%\n', err_C);
if err_C < 1e-4
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

% =========================================================================
% TESTE D — VELOCIDADE TERMINAL COM ARRASTO
% Valida: aero_model.m retorna D = ½ρv²·Sref·CD0 corretamente
%
% Na velocidade terminal: D = m·g → aceleração líquida = 0
% v_t = √(2·m·g / ρ·Sref·CD0)
%
% Duas verificações:
%   Estática: calcula D diretamente e compara com m·g
%   Dinâmica: inicia em v_t e verifica que vz não muda em 5 s
% =========================================================================
fprintf('TESTE D — Velocidade terminal com arrasto\n');
fprintf('  Valida: aero_model computa D = 0.5*rho*v^2*Sref*CD0\n');
fprintf('  Referencia: Anderson (2015) cap. 5; Sutton & Biblarz (2017) cap. 2\n');

% Velocidade terminal teórica (onde arrasto equilibra o peso).
v_t = sqrt(2 * p.m0 * p.g / (p.rho * p.Sref * p.CD0));

fprintf('  v_terminal = sqrt(2*m*g / rho*Sref*CD0) = %.4f m/s\n', v_t);

% Verificação estática
% Avalia arrasto e peso a v_t e confere que a aceleração resultante ~0.
R_vert   = rotation_matrix(0, 0, 0);          % = eye(3)
vel_term = [0; 0; -v_t];                       % descendo em Z_B

F_aero_D = aero_model(vel_term, R_vert, p);
F_grav_D = gravity_body(p.m0, R_vert, p);
F_total  = F_aero_D + F_grav_D;
az_stat  = F_total(3) / p.m0;

D_calc   = norm(F_aero_D);
D_esper  = p.m0 * p.g;
err_D    = abs(D_calc - D_esper) / D_esper * 100;

fprintf('  [Estatico] Arrasto calculado: %.4f N\n', D_calc);
fprintf('  [Estatico] Peso (m*g):        %.4f N\n', D_esper);
fprintf('  [Estatico] Erro D vs mg:      %.8f%%\n', err_D);
fprintf('  [Estatico] az resultante:     %.2e m/s2 (esperado ~0)\n', az_stat);

% Verificação dinâmica
% Inicia exatamente em v_t; se o modelo está certo, vz não muda em 5 s.
x0_D = [0; 0; 5000;   0; 0; -v_t;   p.m0];

[t_D, x_D] = ode45(@(t,x) ode_3dof(t, x, euler_fixo, 0, p), ...
                    [0, 5], x0_D, opts);

vz_var = max(x_D(:,6)) - min(x_D(:,6));

fprintf('  [Dinamico] Variacao de vz em 5s: %.2e m/s (esperado ~0)\n', vz_var);

if err_D < 0.01 && abs(az_stat) < 1e-6 && vz_var < 0.01
    fprintf('  PASSOU\n\n');
else
    fprintf('  FALHOU\n\n');
end

fprintf('=========================================================\n');
fprintf('  RESUMO\n');
fprintf('  Teste B: cinemática 3D — u_B integra em x_I corretamente\n');
fprintf('  Teste C: energia conservada com movimento lateral\n');
fprintf('  Teste D: aero_model computa arrasto corretamente\n');
fprintf('\n');
fprintf('  PROXIMA ETAPA: test_6dof.m\n');
fprintf('  La: rotacao, TVC, acoplamento, CG movel, simetria, torque.\n');
fprintf('=========================================================\n\n');

% ── PLOTS ─────────────────────────────────────────────────────────────────

figure('Name','TESTE B — Projetil Lateral (valida cinematica 3D)', ...
       'NumberTitle','off','Position',[50 50 1000 400]);

% Trajetória parabólica no plano XZ: analítico vs. 3-DOF.
subplot(1,3,1);
plot(x_anal, z_anal, 'r--', 'LineWidth', 2); hold on;
plot(x_B(:,1), x_B(:,3), 'b-', 'LineWidth', 1.5); grid on;
xlabel('x_I [m]'); ylabel('z_I [m]');
legend('Analitico','3-DOF'); axis equal;
title('Teste B — Trajetoria Parabolica (plano XZ)');

% Erro em x_I ao longo do tempo.
subplot(1,3,2);
plot(t_vec, abs(xI_num - x_anal), 'k-', 'LineWidth', 1.5); grid on;
xlabel('Tempo [s]'); ylabel('|Erro x_I| [m]');
title('Teste B — Erro em x_I');

% Conservação de energia 3D (Teste C).
subplot(1,3,3);
plot(t_B, e_num - e_num(1), 'b-', 'LineWidth', 1.5); grid on;
yline(0, 'r--', 'Esperado = 0', 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('\DeltaE [J/kg]');
title('Teste C — Conservacao de Energia 3D');

figure('Name','TESTE D — Velocidade Terminal (valida arrasto)', ...
       'NumberTitle','off','Position',[50 530 700 400]);

% Velocidade vertical estabilizada na terminal.
subplot(1,2,1);
plot(t_D, x_D(:,6), 'b-', 'LineWidth', 2); grid on;
yline(-v_t, 'r--', sprintf('v_t = %.1f m/s', -v_t), 'LineWidth', 1.2);
xlabel('Tempo [s]'); ylabel('w_B [m/s]');
title('Velocidade vertical (deve ser constante)');

% Altitude caindo linearmente (sinal de velocidade constante).
subplot(1,2,2);
plot(t_D, x_D(:,3), 'b-', 'LineWidth', 2); grid on;
xlabel('Tempo [s]'); ylabel('z_I [m]');
title('Altitude (queda linear = terminal confirmado)');

% ── FUNÇÃO LOCAL ──────────────────────────────────────────────────────────

function xdot = ode_3dof(~, x, euler, T, p)
% ODE 3-DOF pura — estado: [x_I; y_I; z_I; u_B; v_B; w_B; m]
% Atitude FIXA (euler é parâmetro, não estado).
% omega = 0 → sem acoplamento rotação-translação.
% Equação de Newton no corpo: v̇_B = F_total_B / m  (sem ω×v_B)

    vel_B = x(4:6);
    m     = x(7);

    % Matriz de rotação da atitude fixa.
    R_IB = rotation_matrix(euler(1), euler(2), euler(3));

    % Forças no corpo: gravidade + arrasto (+ empuxo se T>0).
    F_grav = gravity_body(m, R_IB, p);
    F_aero = aero_model(vel_B, R_IB, p);

    if T > 0
        F_tvc = [0; 0; T];
        mdot  = -T / (p.Isp * p.g0);
    else
        F_tvc = [0; 0; 0];
        mdot  = 0;
    end

    % Cinemática (corpo→inercial) e Newton; massa varia com o empuxo.
    pos_dot = R_IB * vel_B;
    vel_dot = (F_grav + F_aero + F_tvc) / m;

    xdot = [pos_dot; vel_dot; mdot];
end