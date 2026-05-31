function xdot = rocket_6dof(t, x, u, p)
% ROCKET_6DOF  Equações de movimento 6-DOF do foguete, fase terminal de pouso.
%
%   Função de derivadas para o integrador: recebe o estado
%   atual e devolve xdot. Reúne os blocos modulares da planta (TVC, aero,
%   gravidade, CG e inércia variáveis) nas equações de Newton-Euler.

% ── VETOR DE ESTADO (13 componentes) ─────────────────────────────────────
%
%   x( 1: 3) = [x_I; y_I; z_I]    posição no inercial        [m]
%   x( 4: 6) = [u_B; v_B; w_B]    velocidade no corpo        [m/s]
%   x( 7: 9) = [phi; theta; psi]  ângulos de Euler       [rad]
%   x(10:12) = [p_B; q_B; r_B]    taxas angulares no corpo   [rad/s]
%   x(13)    = m                   massa                      [kg]
%
% ── ENTRADAS DE CONTROLE (3 componentes) ─────────────────────────────────
%
%   u(1) = T         empuxo [N]
%   u(2) = delta_y   deflexão TVC em Y [rad]  (>0 → nariz para +Y_B)
%   u(3) = delta_x   deflexão TVC em X [rad]  (>0 → nariz para +X_B)
%
% ── EQUAÇÕES DO MOVIMENTO ─────────────────────────────────────────────────
%
%   EQ1 — Cinemática translacional:
%     ṗ_I = R_IB · v_B
%
%   EQ2 — Dinâmica translacional (Newton no corpo):
%     v̇_B = (1/m) · F_total_B  −  ω_B × v_B
%     O termo (ω_B × v_B) é o acoplamento de Coriolis.
%     F_total_B = F_tvc + F_aero + F_grav   (todos no referencial B)
%
%   EQ3 — Cinemática rotacional (Euler):
%      euler_dot = L(euler) · ω_B
%
%   EQ4 — Dinâmica rotacional (Euler):
%     ω̇_B = J^-1 · (M_tvc  −  ω_B × (J · ω_B))
%     O termo (ω_B × J·ω_B) é o torque giroscópico.
%     Com J diagonal: J^-1 = diag(1/Jxx, 1/Jyy, 1/Jzz)
%
%   EQ5 — Variação de massa (Sutton & Biblarz):
%     m_ponto = −T / (Isp · g0)
% ENTRADAS
%   t — tempo [s]  (não usado explicitamente, mas exigido pelo ode45)
%   x — vetor de estado (13×1)
%   u — entradas de controle (3×1)
%   p — struct de parâmetros (falcon9_params)
%
% SAÍDA
%   xdot — derivadas do estado (13×1)
%
% REFERÊNCIAS: Tewari (2007); Jenie et al. (2019).

% ── 1. DESEMPACOTAR ESTADO ────────────────────────────────────────────────
pos_I   = x(1:3);      % posição inercial
vel_B   = x(4:6);      % velocidade no corpo
euler   = x(7:9);      % ângulos de Euler
omega_B = x(10:12);    % taxas angulares no corpo
m       = x(13);       % massa

phi   = euler(1);
theta = euler(2);
psi   = euler(3);

% ── 2. DESEMPACOTAR E SATURAR ENTRADAS ───────────────────────────────────
T       = u(1);
delta_y = u(2);
delta_z = u(3);

% Saturações físicas (garantia de robustez numérica)
T       = max(p.T_min,     min(p.T_max,     T));
delta_y = max(-p.delta_max, min(p.delta_max, delta_y));
delta_z = max(-p.delta_max, min(p.delta_max, delta_z));

% Proteger massa contra subfluxo
m = max(p.mf, m);

% ── 3. PARÂMETROS DEPENDENTES DA MASSA ───────────────────────────────────
L_TVC          = cg_model(m, p);              % braço de momento TVC [m]
[Jxx, Jyy, Jzz] = inertia_model(m, p);       % tensor de inércia [kg·m2]

% Tensor J e sua inversa (diagonal)
J    = diag([Jxx; Jyy; Jzz]);
Jinv = diag([1/Jxx; 1/Jyy; 1/Jzz]);

% ── 4. MATRIZ DE ROTAÇÃO ──────────────────────────────────────────────────
R_IB = rotation_matrix(phi, theta, psi);  % corpo para inercial

% ── 5. FORÇAS NO REFERENCIAL DO CORPO ────────────────────────────────────
% Todos os vetores de força precisam estar no MESMO referencial (B)
% antes de serem somados — este é o erro mais comum em 6-DOF.

[F_tvc, M_tvc] = tvc_model(T, delta_y, delta_z, L_TVC);
F_aero         = aero_model(vel_B, R_IB, p);
F_grav         = gravity_body(m, R_IB, p);

F_total_B = F_tvc + F_aero + F_grav;   % soma no referencial B ✓

% ── 6. EQ1 — CINEMÁTICA TRANSLACIONAL ────────────────────────────────────
% Derivada da posição inercial = velocidade inercial
% v_I = R_IB * v_B   (rotaciona velocidade do corpo para inercial)
pos_dot = R_IB * vel_B;

% ── 7. EQ2 — DINÂMICA TRANSLACIONAL (Newton no corpo) ─────────────────────
% F = m*a_I, mas reescrito para a velocidade no corpo (Aula 21, Eq. 2):
%   m * v̇_B = F_total_B  −  m * ω_B × v_B
%   v̇_B = (1/m) * F_total_B  −  ω_B × v_B
%
% O termo (-ω_B × v_B) surge da derivada do vetor girante (Eq. 2 da Aula 21).
% Sem ele, o modelo tem dinâmica errada em situações com rotação+translação.
vel_dot = (1/m) * F_total_B - cross(omega_B, vel_B);

% ── 8. EQ3 — CINEMÁTICA ROTACIONAL (Euler) ────────────────────────────────
% Relaciona taxas angulares do corpo com derivadas dos ângulos de Euler.
% Singular em |theta| = 90° — fora do envelope operacional de pouso.
euler_dot = kinematics_euler(euler, omega_B);

% ── 9. EQ4 — DINÂMICA ROTACIONAL (Euler) ──────────────────────────────────
% J * ω̇_B = M_tvc  −  ω_B × (J * ω_B)
%   ω̇_B = J⁻¹ * [M_tvc  −  ω_B × (J * ω_B)]
%
% O termo (-ω_B × J·ω_B) é o torque giroscópico (conservação de momento
% angular no referencial girante). Sem ele, o modelo não reproduz precessão.
% Momento aerodinâmico desprezado (< 0.1% do momento do TVC na fase terminal).
omega_dot = Jinv * (M_tvc - cross(omega_B, J * omega_B));

% ── 10. EQ5 — VARIAÇÃO DE MASSA ─────────────────────────────────────────
% ṁ = -T / (Isp * g₀)   [Sutton & Biblarz, 2017, Eq. 2-3]
if m > p.mf + 1.0    % 1 kg de margem para evitar oscilação
    m_dot = -T / (p.Isp * p.g0);
else
    m_dot = 0.0;     % sem mais propelente
end

% ── 11. MONTAR XDOT ──────────────────────────────────────────────────────
xdot = [pos_dot;    % d/dt [x_I; y_I; z_I]    (3)
        vel_dot;    % d/dt [u_B; v_B; w_B]    (3)
        euler_dot;  % d/dt [phi; theta; psi]  (3)
        omega_dot;  % d/dt [p; q; r]          (3)
        m_dot];     % d/dt [m]                (1)
                    % total: 13 componentes

end
