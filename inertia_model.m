function [Jxx, Jyy, Jzz] = inertia_model(m, p)
% INERTIA_MODEL  Tensor de inércia diagonal em função da massa atual.
%
%   [Jxx, Jyy, Jzz] = inertia_model(m, p)
%
%   O tensor é DIAGONAL (produtos de inércia = 0) pois o foguete é
%   modelado como corpo de revolução axissimétrico em torno de Z_B.
%
%   ATRIBUIÇÃO DOS EIXOS (crítico — nunca confundir):
%     Jzz → rotação em torno de Z_B (eixo de SIMETRIA = ROLL) = menor
%     Jxx → rotação em torno de X_B (eixo LATERAL = PITCH)   = maior
%     Jyy → rotação em torno de Y_B (eixo LATERAL = YAW) = Jxx = maior
%
%   Equações de Euler com tensor diagonal:
%     Jxx * p_dot = Mx - (Jzz - Jyy) * q * r
%     Jyy * q_dot = My - (Jxx - Jzz) * p * r
%     Jzz * r_dot = Mz - (Jyy - Jxx) * p * q
%
%   Modelo: interpolação linear entre m0 e mf.
%
% ENTRADA
%   m — massa atual [kg]
%   p — struct de parâmetros (falcon9_params)
%
% SAÍDA
%   Jxx — momento de inércia em X_B (pitch) [kg·m²]
%   Jyy — momento de inércia em Y_B (yaw)   [kg·m²]  = Jxx
%   Jzz — momento de inércia em Z_B (roll)  [kg·m²]

frac  = (m - p.mf) / (p.m0 - p.mf);
frac  = max(0.0, min(1.0, frac));

Jzz   = p.Jzz_f + (p.Jzz0 - p.Jzz_f) * frac;   % roll
Jxx   = p.Jxx_f + (p.Jxx0 - p.Jxx_f) * frac;   % pitch
Jyy   = Jxx;                                   % yaw = pitch (simetria)

end
