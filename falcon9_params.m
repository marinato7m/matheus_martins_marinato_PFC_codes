function p = falcon9_params()
% =========================================================================
% FALCON9_PARAMS  —  Parâmetros do Falcon 9 Block 5
%
% CONVENCAO DE EIXOS DO CORPO (B) — NUNCA ALTERAR
% -----------------------------------------------------------------------
%   Z_B = eixo de simetria do foguete, aponta para o NARIZ (para cima)
%   X_B = eixo lateral 1  (perpendicular ao eixo de simetria)
%   Y_B = eixo lateral 2  (mao direita com Z_B e X_B)
%
%   Foguete vertical pousando  →  atitude (phi=0, theta=0, psi=0)
%   Empuxo positivo em +Z_B   →  desaceleracao (vz_B negativo = descendo)
%
% CONVENCAO DE EIXOS INERCIAIS (I)
%   Origem = ponto de touchdown
%   Z_I    = para CIMA (altitude positiva)
%   X_I    = Norte,  Y_I = Leste
%
% TENSOR DE INERCIA — HIPOTESE FUNDAMENTAL
% -----------------------------------------------------------------------
%   O foguete e modelado como corpo de revolucao com distribuicao de
%   massa axissimetrica em torno de Z_B. Portanto:
%     - Todos os produtos de inercia sao ZERO (Jxy = Jxz = Jyz = 0)
%     - O tensor J e DIAGONAL no referencial do corpo
%     - X_B, Y_B, Z_B sao os eixos PRINCIPAIS de inercia
%
%   J = diag(Jxx, Jyy, Jzz)
%
%   As equacoes de Euler ficam na forma mais simples:
%     Jxx*p_dot = Mx - (Jzz - Jyy)*q*r    (pitch)
%     Jyy*q_dot = My - (Jxx - Jzz)*p*r    (yaw)
%     Jzz*r_dot = Mz - (Jyy - Jxx)*p*q    (roll)
%
%   Hipotese declarada no TCC: distribuicao axissimetrica; produtos de
%   inercia nulos. Efeitos de assimetria (Octaweb, pernas) desprezados.
%
% FORMULAS DE INERCIA (Meriam & Kraige [R6], Apendice D)
% -----------------------------------------------------------------------
%   Eixo de SIMETRIA (spin / roll):
%     Cilindro OCO (casca fina):   I_sim = m * R^2
%     Cilindro SOLIDO:             I_sim = (1/2) * m * R^2
%
%   Eixo LATERAL (tombar / pitch / yaw):
%     Cilindro OCO (casca fina):   I_lat = m * (R^2/2  + L^2/12)
%     Cilindro SOLIDO:             I_lat = m * (R^2/4  + L^2/12)
%
%   Em nossa convencao:
%     Jzz (roll)      ←→  I_sim  (eixo Z_B = simetria)  → PEQUENO
%     Jxx = Jyy       ←→  I_lat  (eixos X_B, Y_B)       → GRANDE
%
% REFERENCIAS
%   [R1] SPACEX. Falcon Payload User's Guide. mai. 2025.
%        <https://www.spacex.com/assets/media/falcon-users-guide-2025-05-09.pdf>
%   [R2] SUTTON, G. P.; BIBLARZ, O. Rocket Propulsion Elements. 9. ed.
%        Wiley, 2017. ISBN 978-1-118-75365-1.
%   [R3] JENIE, Y. I.; SUARJAYA, W. W. H.; POETRO, R. E. Falcon 9 Rocket
%        Launch Modeling and Simulation with Thrust Vectoring Control and
%        Scheduling. IEEE ACDT 2019. DOI: 10.1109/ACDT47198.2019.9072837.
%   [R4] TEWARI, A. Atmospheric and Space Flight Dynamics. Birkhauser, 2007.
%   [R5] ANDERSON, J. D. Introduction to Flight. 8. ed. McGraw-Hill, 2015.
%   [R6] MERIAM, J. L.; KRAIGE, L. G. Engineering Mechanics: Dynamics.
%        8. ed. Wiley, 2015.
% =========================================================================

% -------------------------------------------------------------------------
%  1. GEOMETRIA
% -------------------------------------------------------------------------
p.L         = 41.2;          % [m]   Comprimento 1a fase  [R3]
p.D         = 3.66;          % [m]   Diametro              [R1]
p.R         = p.D / 2;       % [m]   Raio
p.Sref      = pi * p.R^2;    % [m2]  Area de referencia aerodinamica

% Comprimento da secao de tanques (para inercia do propelente)
% Tanques LOX+RP-1 ocupam ~75% do comprimento da 1a fase.
p.L_tank    = 0.75 * p.L;    % [m]   ~30.9m  [C]

% -------------------------------------------------------------------------
%  2. MASSAS
% -------------------------------------------------------------------------
p.m_seco    = 22200;         % [kg]  Massa seca (estrutura + motores)  [R3]
p.m_motores = 9 * 490;       % [kg]  9 motores Merlin 1D, ~490 kg cada [R2]
p.m0        = 30000;         % [kg]  Massa ao inicio do landing burn   [R3] 30 000
p.m_prop0   = p.m0 - p.m_seco; % [kg]  Propelente inicial = 7800 kg
p.mf        = 22500;         % [kg]  Massa no touchdown                [C]

assert(p.m0  > p.mf,       'ERRO: massa final maior que inicial');
assert(p.mf >= p.m_seco,   'ERRO: massa final menor que massa seca');

% -------------------------------------------------------------------------
%  3. PROPULSAO (1 motor Merlin 1D Block 5, nivel do mar)
% -------------------------------------------------------------------------
p.T_max        = 845000;     % [N]  Empuxo maximo SL            [R1]
p.throttle_min = 0.40;       % [-]  Throttle minimo             [R1][R3]
p.T_min        = p.throttle_min * p.T_max;   % = 338000 N
p.Isp          = 282;        % [s]  Isp SL (Merlin 1D)          [R2]
p.g0           = 9.80665;    % [m/s2] ISO 80000-3
p.mdot_max     = p.T_max / (p.Isp * p.g0);  % ~305.6 kg/s
p.mdot_min     = p.T_min / (p.Isp * p.g0);  % ~122.2 kg/s
p.n_mot_pouso  = 1;          % [-]  Landing burn: 1 motor central [R1]

% -------------------------------------------------------------------------
%  4. TVC
% -------------------------------------------------------------------------
p.delta_max     = deg2rad(6.0);   % [rad]   Deflexao max. gimbal  [R3]
p.delta_dot_max = deg2rad(20.0);  % [rad/s] Taxa max. deflexao    [C][R4]

% -------------------------------------------------------------------------
%  5. POSICAO DO CG (medida da BASE — bocal dos motores)
% -------------------------------------------------------------------------
%
% O CG DESCE conforme o propelente e consumido.
%
% Por que desce?
%   - Motores Merlin (4410 kg) ficam FIXOS na base (~0.5m da base)
%   - Propelente (que queima) fica nos tanques ACIMA (~24.7m da base)
%   - Ao consumir propelente, o CG migra em direcao aos motores (base)
%
% Calculo por componentes:
%   Em m0 = 30000 kg:
%     r_CG = (4410*0.5 + 17790*20.6 + 7800*24.7) / 30000 ~ 18.5 m
%   Em mf = 22500 kg (propelente ~ 300 kg):
%     r_CG = (4410*0.5 + 17790*20.6 + 300*24.7)  / 22500 ~ 17.0 m
%
% Impacto no controle:
%   L_TVC = r_CG  (distancia CG-bocal = braco do momento do TVC)
%   L_TVC diminui ~1.5 m → autoridade de controle reduz ao final do pouso

p.r_CG_0 = 18.5;  % [m]  CG em m = m0 (inicio do landing burn)  [C]
p.r_CG_f = 17.0;  % [m]  CG em m = mf (touchdown)               [C]

% Modelo linear: ver funcao cg_model() ao final do arquivo.

% -------------------------------------------------------------------------
%  6. TENSOR DE INERCIA  (DIAGONAL — produtos de inercia nulos)
% -------------------------------------------------------------------------
% --- Inercia da estrutura (constante, pois m_seco nao varia) ---
Jzz_casca = p.m_seco * p.R^2;
Jxx_casca = p.m_seco * (p.R^2/2 + p.L^2/12);

% --- Inercia do propelente em m = m0 ---
Jzz_prop0 = 0.5 * p.m_prop0 * p.R^2;
Jxx_prop0 = p.m_prop0 * (p.R^2/4 + p.L_tank^2/12);

% --- Inercia total em m = m0 ---
p.Jzz0 = Jzz_casca + Jzz_prop0;   % roll
p.Jxx0 = Jxx_casca + Jxx_prop0;   % pitch
p.Jyy0 = p.Jxx0;                  % yaw = pitch (simetria)

% --- Inercia do propelente em m = mf ---
m_prop_f  = p.mf - p.m_seco;      % = 300 kg
Jzz_prop_f = 0.5 * m_prop_f * p.R^2;
Jxx_prop_f = m_prop_f * (p.R^2/4 + p.L_tank^2/12);

% --- Inercia total em m = mf ---
p.Jzz_f = Jzz_casca + Jzz_prop_f;
p.Jxx_f = Jxx_casca + Jxx_prop_f;
p.Jyy_f = p.Jxx_f;

% Modelo linear: ver funcao inertia_model() ao final do arquivo.

% -------------------------------------------------------------------------
%  7. AERODINAMICA
% -------------------------------------------------------------------------
% Fase terminal: v < 70 m/s, M < 0.21, angulos de ataque < 5 deg
% → sustentacao nula, apenas arrasto  [R4]
p.CD0      = 0.30;           % [-]     Coef. arrasto subsonico  [R4][R5]
p.CL_alpha = 2.0;            % [1/rad] Teoria corpo esbelto     [R5]
p.rho      = 1.225;          % [kg/m3] ISA, nivel do mar (ISO 2533)
p.a_som    = 340.3;          % [m/s]   ISA, T = 288.15 K

% -------------------------------------------------------------------------
%  8. GRAVIDADE
% -------------------------------------------------------------------------
p.g   = 9.80665;             % [m/s2]   ISO 80000-3
p.g_I = [0; 0; -p.g];       % [m/s2]   vetor no inercial (Z_I para cima)

% -------------------------------------------------------------------------
%  9. CONDICOES INICIAIS
% -------------------------------------------------------------------------
% --- Canal vertical (z) --------------------------------------------------
p.h0     = 850;              % [m]    Altitude inicial landing burn  [C]
p.vz0    = -50.0;            % [m/s]  Velocidade vertical (neg=descendo)
p.az0    = 0.0;              % [m/s2]
p.hf     = 0.0;              % [m]    Altitude de touchdown
p.vzf    = -2.0;              % [m/s]  Velocidade final desejada
p.azf    = 4.0;              % [m/s2]

% --- Canal lateral x -----------------------------------------------------
p.x0     = 20.0;             % [m]    Offset lateral inicial em x 
p.vx0    = -3.0;              % [m/s]  Velocidade lateral inicial em x 
p.ax0    = 0.0;              % [m/s2] Aceleracao lateral inicial em x
p.xf     = 0.0;              % [m]    Alvo final em x (pouso centrado)
p.vxf    = 0.0;              % [m/s]  Velocidade final em x
p.axf    = 0.0;              % [m/s2] Aceleracao final em x

% --- Canal lateral y -----------------------------------------------------
p.y0     = -8.0;            % [m]    Offset lateral inicial em y
p.vy0    = 5.0;              % [m/s]  Velocidade lateral inicial em y
p.ay0    = 0.0;              % [m/s2] Aceleracao lateral inicial em y
p.yf     = 0.0;              % [m]    Alvo final em y
p.vyf    = 0.0;              % [m/s]  Velocidade final em y
p.ayf    = 0.0;              % [m/s2] Aceleracao final em y

% --- Atitude e taxas angulares -------------------------------------------
p.phi0   = deg2rad(-5.0);              % [rad]  Pitch inicial
p.theta0 = deg2rad(2.0);     % [rad]  Yaw inicial
p.psi0   = 0.0;              % [rad]  Roll inicial
p.p0     = 0.0;              % [rad/s] Taxa de Pitch
p.q0     = deg2rad(0.5);     % [rad/s] Taxa de Yaw
p.r0     = 0.0;              % [rad/s] Taxa de Roll

% --- Ambiente ------------------------------------------------------------
p.v_vento = [0; 0; 0];       % [m/s]  Vento nominal zero; variar no MC
% -------------------------------------------------------------------------
%  10. SOLVER
% -------------------------------------------------------------------------
p.dt          = 0.001;       % [s]   Passo fixo
p.t_max       = 30.0;        % [s]   Tempo maximo de simulacao
p.h_touchdown = 0.10;        % [m]   Altitude de deteccao de touchdown
p.v_touchdown = 1.0;         % [m/s] Velocidade maxima aceitavel em TD

% -------------------------------------------------------------------------
%  11. LIMITES DE DESEMPENHO
% -------------------------------------------------------------------------
p.az_max = (p.T_max - p.m0*p.g) / p.m0;  % ~ +18.4 m/s2
p.az_min = (p.T_min - p.m0*p.g) / p.m0;  % ~  +1.5 m/s2
% =========================================================================
%  RESUMO IMPRESSO
% =========================================================================
fprintf('\n=== PARAMETROS FALCON 9 BLOCK 5 — FASE TERMINAL ===\n');
fprintf('Geometria:\n');
fprintf('  L (1a fase) : %.1f m   D : %.2f m   Sref : %.2f m2\n', ...
        p.L, p.D, p.Sref);
fprintf('CG (da base):\n');
fprintf('  em m0 = %.0f kg : r_CG = %.1f m (%.0f%% de L)\n', ...
        p.m0, p.r_CG_0, 100*p.r_CG_0/p.L);
fprintf('  em mf = %.0f kg : r_CG = %.1f m (%.0f%% de L)\n', ...
        p.mf, p.r_CG_f, 100*p.r_CG_f/p.L);
fprintf('  CG desce %.1f m durante o landing burn\n', p.r_CG_0 - p.r_CG_f);
fprintf('Massas:\n');
fprintf('  m_seco = %.0f kg   m0 = %.0f kg   m_prop0 = %.0f kg   mf = %.0f kg\n',...
        p.m_seco, p.m0, p.m_prop0, p.mf);
fprintf('Inercia em m0 (casca oca + propelente solido):\n');
fprintf('  Jzz (ROLL  — Z_B = simetria) : %.3e kg.m2  [PEQUENO]\n', p.Jzz0);
fprintf('  Jxx (PITCH — X_B = lateral)  : %.3e kg.m2  [GRANDE]\n',  p.Jxx0);
fprintf('  Jyy (YAW   — Y_B = lateral)  : %.3e kg.m2  [GRANDE]\n',  p.Jyy0);
fprintf('  Razao Jxx/Jzz : %.0f  (corpo esbelto — correto)\n', p.Jxx0/p.Jzz0);
fprintf('Inercia em mf:\n');
fprintf('  Jzz_f : %.3e   Jxx_f : %.3e   Jyy_f : %.3e  kg.m2\n',...
        p.Jzz_f, p.Jxx_f, p.Jyy_f);
fprintf('  Obs: todos os produtos de inercia = 0 (corpo de revolucao)\n');
fprintf('Propulsao:\n');
fprintf('  T_max = %.0f kN   T_min = %.0f kN   Isp = %.0f s\n',...
        p.T_max/1e3, p.T_min/1e3, p.Isp);
fprintf('  mdot_max = %.1f kg/s\n', p.mdot_max);
fprintf('TVC:\n');
fprintf('  delta_max = %.1f deg   dot_max = %.1f deg/s\n',...
        rad2deg(p.delta_max), rad2deg(p.delta_dot_max));
fprintf('Desempenho (m = m0):\n');
fprintf('  az_max = %+.2f m/s2   az_min = %+.2f m/s2\n', p.az_max, p.az_min);
fprintf('===================================================\n\n');



if evalin('base','exist(''v_vento_mc'',''var'')')
    p.v_vento = evalin('base','v_vento_mc');
end
end


% =========================================================================
%  FUNCOES AUXILIARES
% =========================================================================

function r_cg = cg_model(m, p)
% CG_MODEL  Posicao do CG em relacao a BASE do foguete (bocal).
%
%   Modelo: interpolacao linear entre r_CG_0 (em m0) e r_CG_f (em mf).
%   O CG desce conforme o propelente queima (motores pesados na base).
%
%   Esta funcao retorna o braco de alavanca do TVC:
%     L_TVC(t) = cg_model(m(t), p)
%
%   As equacoes de Newton/Euler sao validas pois sao escritas para o
%   CG instantaneo — mover o CG dentro do corpo nao viola nenhuma lei.
%
% ENTRADA: m [kg], p [struct]
% SAIDA:   r_cg [m] — distancia CG-bocal

    frac  = (m - p.mf) / (p.m0 - p.mf);
    frac  = max(0, min(1, frac));     % saturacao
    r_cg  = p.r_CG_f + (p.r_CG_0 - p.r_CG_f) * frac;
end


function [Jxx, Jyy, Jzz] = inertia_model(m, p)
% INERTIA_MODEL  Tensor de inercia diagonal em funcao da massa atual.

% ENTRADA: m [kg], p [struct]
% SAIDA:   Jxx, Jyy, Jzz [kg.m2]

    frac  = (m - p.mf) / (p.m0 - p.mf);
    frac  = max(0, min(1, frac));
    Jzz   = p.Jzz_f + (p.Jzz0 - p.Jzz_f) * frac;
    Jxx   = p.Jxx_f + (p.Jxx0 - p.Jxx_f) * frac;
    Jyy   = Jxx;                                
end

% =========================================================================
%  TABELA CONSOLIDADA
% =========================================================================
%
%  Parametro      Valor            Confianca  Fonte
%  ------------------------------------------------------------------
%  L              41.2 m           [A]        [R3] Jenie et al. 2019
%  D              3.66 m           [A]        [R1] SpaceX User's Guide 2025
%  m_seco         22.200 kg        [A]        [R3] Jenie et al. 2019
%  m0             30.000 kg        [A]        [R3] Jenie et al. 2019
%  T_max          845 kN           [A]        [R1] SpaceX User's Guide 2025
%  Throttle min.  40%              [A]        [R1]+[R3]
%  Isp_SL         282 s            [A]        [R2] Sutton & Biblarz, 9a ed.
%  delta_max      +/-6 deg         [A]        [R3] Jenie et al. 2019
%  r_CG_0         18.5 m           [C]        Calculo por componentes
%  r_CG_f         17.0 m           [C]        Calculo por componentes
%  Jyy0           ver impressao    [C]        [R6] Meriam & Kraige, formula
%  CD0            0.30             [B]        [R4][R5]
%  CL_alpha       2.0 /rad         [B]        [R5] Anderson, cap. 5
%  rho            1.225 kg/m3      [A]        ISO 2533
%  g              9.80665 m/s2     [A]        ISO 80000-3
%
%  Parametros NAO publicados pela SpaceX (hipoteses no TCC):
%    r_CG, J_yy  -> modelo geometrico; 
% =========================================================================