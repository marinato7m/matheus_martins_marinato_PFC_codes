function traj = gerar_traj_rapido(h0, vz0, x0, y0, u_B0, v_B0, phi0, theta0, p)
% GERAR_TRAJ_RAPIDO  Geração rápida de trajetória S-curve 3D.
%
%
%   Replica a lógica do run_scurve, mas substitui o drop test Simulink
%   por queda livre analítica. Isso reduz o tempo de ~5s para <0.1s,
%   viabilizando o uso dentro de loops de otimização.
%
%   ENTRADAS:
%     h0, vz0       — altitude [m] e vel. vertical [m/s] (body, negativa = descendo)
%     x0, y0        — posição lateral inercial [m]
%     u_B0, v_B0    — velocidade lateral no body frame [m/s]
%     phi0, theta0  — atitude inicial [rad]
%     p             — struct de falcon9_params()
%
%   SAÍDA:
%     traj — struct idêntico ao do run_scurve (com timeseries para Simulink)
%            Retorna [] se a manobra é impossível para estas CIs.
%
%   SIMPLIFICAÇÕES vs run_scurve v6:
%     - Queda livre analítica (sem arrasto aerodinâmico)
%     - Atitude constante durante queda livre
%     - alpha_ign = 0.5 (meio da janela)
%
% Matheus Martins Marinato — TCC UFSM 2025/2026

    g            = p.g;
    dt           = p.dt;
    m0           = p.m0;
    empuxo_max   = p.T_max;
    T_min_motor  = p.T_max * p.throttle_min;
    margem       = 0.90;
    
    hf  = p.hf;    vzf = p.vzf;   azf = p.azf;
    xf  = 0;       vxf = 0;       axf = 0;
    yf  = 0;       vyf = 0;       ayf = 0;
    
    psi0 = 0;
    
    % ── Converter velocidades body para inercial ───────────────────────────
    R_IB  = rotation_matrix(phi0, theta0, psi0);
    V_I   = R_IB * [u_B0; v_B0; vz0];
    vx0_I = V_I(1);
    vy0_I = V_I(2);
    vz0_I = V_I(3);
    
    % ══════════════════════════════════════════════════════════════════════
    %  FASE 0 — QUEDA LIVRE ANALÍTICA
    % ══════════════════════════════════════════════════════════════════════
    
    % Tempo de impacto (queda livre, h(t)=0)
    % h0 + vz0_I*t - 0.5*g*t² = 0  ->  t = (-vz0_I + sqrt(vz0_I² + 2*g*h0)) / g
    discriminante = vz0_I^2 + 2*g*h0;
    if discriminante < 0
        traj = [];
        return;
    end
    t_impacto = (-vz0_I + sqrt(discriminante)) / g;
    
    % Vetor temporal da queda (resolução 0.01s como o Simulink)
    t_drop = (0 : 0.01 : t_impacto)';
    
    % Estados analíticos durante queda livre
    drop_h   = h0    + vz0_I * t_drop - 0.5 * g * t_drop.^2;
    drop_vz  = vz0_I - g * t_drop;
    drop_x   = x0    + vx0_I * t_drop;
    drop_y   = y0    + vy0_I * t_drop;
    drop_vx  = vx0_I * ones(size(t_drop));
    drop_vy  = vy0_I * ones(size(t_drop));
    drop_phi   = phi0   * ones(size(t_drop));   % atitude constante
    drop_theta = theta0 * ones(size(t_drop));
    
    % ══════════════════════════════════════════════════════════════════════
    %  FASE 1 — SCAN RÁPIDO DA JANELA DE IGNIÇÃO
    % ══════════════════════════════════════════════════════════════════════
    
    SCAN_DT  = 0.2;          % resolução do scan [s]
    T_VIOL_TOL = 0.03;       % 3% tolerância
    
    scan_step = max(1, round(SCAN_DT / 0.01));
    scan_idx  = 1 : scan_step : length(t_drop);
    N_scan    = length(scan_idx);
    
    scan_t      = NaN(N_scan, 1);
    scan_h_min  = NaN(N_scan, 1);
    scan_T_viol = NaN(N_scan, 1);
    scan_valid  = false(N_scan, 1);
    
    for k = 1:N_scan
        ii = scan_idx(k);
        
        t_i   = t_drop(ii);
        h_i   = drop_h(ii);
        x_i   = drop_x(ii);
        y_i   = drop_y(ii);
        vx_i  = drop_vx(ii);
        vy_i  = drop_vy(ii);
        vz_i  = drop_vz(ii);
        phi_i = drop_phi(ii);
        th_i  = drop_theta(ii);
        
        scan_t(k) = t_i;
        
        if h_i <= hf + 0.5
            continue;
        end
        
        % Aceleração na ignição (T_min na atitude atual)
        ax_i = (T_min_motor / m0) * sin(th_i)  * cos(phi_i);
        ay_i = -(T_min_motor / m0) * sin(phi_i);
        az_i = (T_min_motor / m0) * cos(th_i)  * cos(phi_i) - g;
        
        [ok_k, h_min_k, T_viol_k] = scan_frame_rapido( ...
            h_i,  vz_i, az_i,  hf, vzf, azf, ...
            x_i,  vx_i, ax_i,  xf, vxf, axf, ...
            y_i,  vy_i, ay_i,  yf, vyf, ayf, ...
            empuxo_max, T_min_motor, m0, margem, dt, g);
        
        scan_valid(k)   = ok_k;
        scan_h_min(k)   = h_min_k;
        scan_T_viol(k)  = T_viol_k;
    end
    
    % ── Determinar janela ───────────────────────────────────────────────
    h_floor = -0.5;
    
    idx_late = find(scan_valid & scan_h_min >= h_floor);
    if isempty(idx_late)
        traj = [];
        return;
    end
    i_late = idx_late(end);
    t_late = scan_t(i_late);
    
    idx_early = find(scan_valid & scan_T_viol <= T_VIOL_TOL);
    if ~isempty(idx_early)
        i_early = idx_early(1);
        t_early = scan_t(i_early);
    else
        [~, ib] = min(scan_T_viol);
        i_early = ib;
        t_early = scan_t(i_early);
    end
    
    if t_late < t_early
        traj = [];
        return;
    end
    
    % ── Selecionar instante de ignição (alpha = 0.5) ───────────────────
    alpha_ign = 0.5;
    t_ign = t_late - alpha_ign * (t_late - t_early);
    
    [~, idx_ign] = min(abs(t_drop - t_ign));
    t_ign = t_drop(idx_ign);
    
    % ══════════════════════════════════════════════════════════════════════
    %  FASE 2 — GERAR S-CURVE DO BURN
    % ══════════════════════════════════════════════════════════════════════
    
    h_ign     = drop_h(idx_ign);
    x_ign     = drop_x(idx_ign);
    y_ign     = drop_y(idx_ign);
    vx_ign    = drop_vx(idx_ign);
    vy_ign    = drop_vy(idx_ign);
    vz_ign    = drop_vz(idx_ign);
    phi_ign   = drop_phi(idx_ign);
    theta_ign = drop_theta(idx_ign);
    
    ax_ign = (T_min_motor / m0) * sin(theta_ign) * cos(phi_ign);
    ay_ign = -(T_min_motor / m0) * sin(phi_ign);
    az_ign = (T_min_motor / m0) * cos(theta_ign) * cos(phi_ign) - g;
    
    try
        traj_burn = scurve_calc_3d( ...
            h_ign, vz_ign, az_ign,  hf, vzf, azf, ...
            x_ign, vx_ign, ax_ign,  xf, vxf, axf, ...
            y_ign, vy_ign, ay_ign,  yf, vyf, ayf, ...
            empuxo_max, T_min_motor, m0, margem, dt, g);
    catch
        traj = [];
        return;
    end
    
    T_burn = traj_burn.T;
    
    % ══════════════════════════════════════════════════════════════════════
    %  FASE 3 — CONCATENAR QUEDA LIVRE + BURN
    % ══════════════════════════════════════════════════════════════════════
    
    % Queda livre até ignição (sem último ponto — duplicado)
    i_ff = 1:(idx_ign - 1);
    
    t_burn_shifted = traj_burn.time + t_ign;
    
    t_total   = [t_drop(i_ff);     t_burn_shifted];
    h_total   = [drop_h(i_ff);     traj_burn.h_d];
    vz_total  = [drop_vz(i_ff);    traj_burn.vz_d];
    az_total  = [-g*ones(length(i_ff),1);  traj_burn.az_d];
    
    x_total   = [drop_x(i_ff);     traj_burn.x_d];
    vx_total  = [drop_vx(i_ff);    traj_burn.vx_d];
    ax_total  = [zeros(length(i_ff),1);    traj_burn.ax_d];
    
    y_total   = [drop_y(i_ff);     traj_burn.y_d];
    vy_total  = [drop_vy(i_ff);    traj_burn.vy_d];
    ay_total  = [zeros(length(i_ff),1);    traj_burn.ay_d];
    
    theta_total = [drop_theta(i_ff);  traj_burn.theta_d];
    phi_total   = [drop_phi(i_ff);    traj_burn.phi_d];
    Tff_total   = [zeros(length(i_ff),1);  traj_burn.Tff_d];
    
    % ── Converter velocidades para body frame ──────────────────────────
    psi_total = zeros(size(t_total));
    u_total   = zeros(size(t_total));
    v_total   = zeros(size(t_total));
    w_total   = zeros(size(t_total));
    
    for i = 1:length(t_total)
        R_IB_i = rotation_matrix(phi_total(i), theta_total(i), psi_total(i));
        R_BI_i = R_IB_i';
        V_B    = R_BI_i * [vx_total(i); vy_total(i); vz_total(i)];
        u_total(i) = V_B(1);
        v_total(i) = V_B(2);
        w_total(i) = V_B(3);
    end
    
    % ══════════════════════════════════════════════════════════════════════
    %  FASE 4 — MONTAR STRUCT TRAJ (idêntico ao run_scurve)
    % ══════════════════════════════════════════════════════════════════════
    
    traj.T       = t_ign + T_burn;
    traj.T_burn  = T_burn;
    traj.t_ign   = t_ign;
    traj.h_ign   = h_ign;
    traj.v_ign   = vz_ign;
    traj.dt      = dt;
    traj.time    = t_total;
    
    traj.coef    = traj_burn.coef;
    traj.s       = traj_burn.s;
    traj.h_d     = h_total;
    traj.vz_d    = vz_total;
    traj.az_d    = az_total;
    traj.h_ref   = timeseries(h_total,  t_total, 'Name', 'h_ref');
    traj.vz_ref  = timeseries(vz_total, t_total, 'Name', 'vz_ref');
    traj.az_ref  = timeseries(az_total, t_total, 'Name', 'az_ref');
    
    traj.coef_x  = traj_burn.coef_x;
    traj.x_d     = x_total;
    traj.vx_d    = vx_total;
    traj.ax_d    = ax_total;
    traj.x_ref   = timeseries(x_total,  t_total, 'Name', 'x_ref');
    traj.vx_ref  = timeseries(vx_total, t_total, 'Name', 'vx_ref');
    traj.ax_ref  = timeseries(ax_total, t_total, 'Name', 'ax_ref');
    
    traj.coef_y  = traj_burn.coef_y;
    traj.y_d     = y_total;
    traj.vy_d    = vy_total;
    traj.ay_d    = ay_total;
    traj.y_ref   = timeseries(y_total,  t_total, 'Name', 'y_ref');
    traj.vy_ref  = timeseries(vy_total, t_total, 'Name', 'vy_ref');
    traj.ay_ref  = timeseries(ay_total, t_total, 'Name', 'ay_ref');
    
    traj.theta_d   = theta_total;
    traj.phi_d     = phi_total;
    traj.Tff_d     = Tff_total;
    traj.theta_ref = timeseries(theta_total, t_total, 'Name', 'theta_ref');
    traj.phi_ref   = timeseries(phi_total,   t_total, 'Name', 'phi_ref');
    traj.Tff_ref   = timeseries(Tff_total,   t_total, 'Name', 'Tff_ref');
    
    traj.u_d    = u_total;
    traj.v_d    = v_total;
    traj.w_d    = w_total;
    traj.u_ref  = timeseries(u_total, t_total, 'Name', 'u_ref');
    traj.v_ref  = timeseries(v_total, t_total, 'Name', 'v_ref');
    traj.w_ref  = timeseries(w_total, t_total, 'Name', 'w_ref');
end


% ═════════════════════════════════════════════════════════════════════════
%  FUNÇÕES LOCAIS (extraídas do run_scurve.m)
% ═════════════════════════════════════════════════════════════════════════

function [ok, h_min, T_viol] = scan_frame_rapido( ...
    h0,vz0,az0, hf,vzf,azf, ...
    x0,vx0,ax0, xf,vxf,axf, ...
    y0,vy0,ay0, yf,vyf,ayf, ...
    empuxo_max, T_min_motor, m0, margem, dt, g)

    ok     = false;
    h_min  = NaN;
    T_viol = NaN;

    dh = h0 - hf;
    if dh <= 0;  return;  end

    try
        f_res = @(T_) max_thrust_dado_T(T_, ...
            h0,vz0,az0, hf,vzf,azf, dh, dt, ...
            x0,vx0,ax0, xf,vxf,axf, ...
            y0,vy0,ay0, yf,vyf,ayf, ...
            m0, g) - (margem * empuxo_max);

        opts = optimset('TolX',1e-3,'Display','off');
        T_ot = fzero(f_res, [5, 300], opts);
    catch
        return;
    end

    try
        [~, t, ~, ~, ~, ~, h_d, ~, az_d] = ...
            resolver_scurve_z(T_ot, h0,vz0,az0,hf,vzf,azf,dh,dt);
    catch
        return;
    end

    h_min = min(h_d);
    ok    = true;

    if h_min < 0
        T_viol = 1.0;
        return;
    end

    try
        [~, ~, ~, ax_d] = resolver_canal(T_ot, x0,vx0,ax0, xf,vxf,axf, t);
        [~, ~, ~, ay_d] = resolver_canal(T_ot, y0,vy0,ay0, yf,vyf,ayf, t);
    catch
        T_viol = 1.0;
        return;
    end

    Tff_d = m0 * sqrt(ax_d.^2 + ay_d.^2 + (g + az_d).^2);
    T_viol = sum(Tff_d < T_min_motor) / length(Tff_d);
end

function traj_burn = scurve_calc_3d( ...
    h0,vz0,az0, hf,vzf,azf, ...
    x0,vx0,ax0, xf,vxf,axf, ...
    y0,vy0,ay0, yf,vyf,ayf, ...
    empuxo_max, T_min_motor, m0, margem, dt, g) %#ok<INUSD>

    dh = h0 - hf;

    f_res = @(T_) max_thrust_dado_T(T_, ...
        h0,vz0,az0, hf,vzf,azf, dh, dt, ...
        x0,vx0,ax0, xf,vxf,axf, ...
        y0,vy0,ay0, yf,vyf,ayf, ...
        m0, g) - (margem * empuxo_max);

    opts = optimset('TolX',1e-4,'Display','off');
    T_ot = fzero(f_res, [5, 300], opts);

    [coef_z_norm, t, ~, s, ~, ~, h_d, vz_d, az_d] = ...
        resolver_scurve_z(T_ot, h0,vz0,az0,hf,vzf,azf,dh,dt);

    [coef_x, x_d, vx_d, ax_d] = resolver_canal( ...
        T_ot, x0, vx0, ax0, xf, vxf, axf, t);

    [coef_y, y_d, vy_d, ay_d] = resolver_canal( ...
        T_ot, y0, vy0, ay0, yf, vyf, ayf, t);

    theta_d = atan2(ax_d, g + az_d);
    phi_d   = -atan2(ay_d, g + az_d);
    Tff_d   = m0 * sqrt(ax_d.^2 + ay_d.^2 + (g + az_d).^2);

    traj_burn.T       = T_ot;
    traj_burn.dt      = dt;
    traj_burn.time    = t;
    traj_burn.coef    = coef_z_norm;
    traj_burn.s       = s;
    traj_burn.h_d     = h_d;      traj_burn.vz_d  = vz_d;    traj_burn.az_d  = az_d;
    traj_burn.coef_x  = coef_x;
    traj_burn.x_d     = x_d;      traj_burn.vx_d  = vx_d;    traj_burn.ax_d  = ax_d;
    traj_burn.coef_y  = coef_y;
    traj_burn.y_d     = y_d;      traj_burn.vy_d  = vy_d;    traj_burn.ay_d  = ay_d;
    traj_burn.theta_d = theta_d;  traj_burn.phi_d = phi_d;   traj_burn.Tff_d = Tff_d;
end

function T_pk = max_thrust_dado_T(T_, ...
    h0,vz0,az0, hf,vzf,azf, dh, dt, ...
    x0,vx0,ax0, xf,vxf,axf, ...
    y0,vy0,ay0, yf,vyf,ayf, ...
    m0, g)

    try
        [~, t, ~, ~, ~, ~, ~, ~, az_d] = ...
            resolver_scurve_z(T_, h0,vz0,az0,hf,vzf,azf,dh,dt);
        [~, ~, ~, ax_d] = resolver_canal(T_, x0,vx0,ax0, xf,vxf,axf, t);
        [~, ~, ~, ay_d] = resolver_canal(T_, y0,vy0,ay0, yf,vyf,ayf, t);
        Tff_d = m0 * sqrt(ax_d.^2 + ay_d.^2 + (g + az_d).^2);
        T_pk = max(Tff_d);
    catch
        T_pk = NaN;
    end
end

function [coef,t,tau,s,s_d1,s_d2,h_d,vz_d,az_d] = ...
         resolver_scurve_z(T_, h0,vz0,az0,hf,vzf,azf,dh,dt)

    s0   =  0;
    s1   =  1;
    sd0  = -vz0 * T_  / dh;
    sd1  = -vzf * T_  / dh;
    sdd0 = -az0 * T_^2 / dh;
    sdd1 = -azf * T_^2 / dh;

    M = [1  0  0   0   0    0  ;
         0  1  0   0   0    0  ;
         0  0  2   0   0    0  ;
         1  1  1   1   1    1  ;
         0  1  2   3   4    5  ;
         0  0  2   6   12   20 ];

    b    = [s0; sd0; sdd0; s1; sd1; sdd1];
    coef = M \ b;

    t   = (0:dt:T_)';
    tau = t / T_;

    s    = coef(1) + coef(2)*tau      + coef(3)*tau.^2   + ...
           coef(4)*tau.^3 + coef(5)*tau.^4 + coef(6)*tau.^5;
    s_d1 =           coef(2)          + 2*coef(3)*tau    + ...
           3*coef(4)*tau.^2 + 4*coef(5)*tau.^3 + 5*coef(6)*tau.^4;
    s_d2 =                     2*coef(3)                 + ...
           6*coef(4)*tau    + 12*coef(5)*tau.^2 + 20*coef(6)*tau.^3;

    h_d  =  h0 - dh * s;
    vz_d = -(dh / T_)   * s_d1;
    az_d = -(dh / T_^2) * s_d2;
end

function [coef, p_d, v_d, a_d] = resolver_canal(T_, p0, v0, a0, pf, vf, af, t)

    tau = t / T_;

    M = [1  0  0   0   0    0  ;
         0  1  0   0   0    0  ;
         0  0  2   0   0    0  ;
         1  1  1   1   1    1  ;
         0  1  2   3   4    5  ;
         0  0  2   6   12   20 ];

    b    = [p0; T_*v0; T_^2*a0; pf; T_*vf; T_^2*af];
    coef = M \ b;

    p_d = coef(1) + coef(2)*tau      + coef(3)*tau.^2   + ...
          coef(4)*tau.^3 + coef(5)*tau.^4 + coef(6)*tau.^5;
    v_d = (coef(2) + 2*coef(3)*tau   + 3*coef(4)*tau.^2 + ...
           4*coef(5)*tau.^3 + 5*coef(6)*tau.^4) / T_;
    a_d = (2*coef(3) + 6*coef(4)*tau + 12*coef(5)*tau.^2 + ...
           20*coef(6)*tau.^3) / T_^2;
end