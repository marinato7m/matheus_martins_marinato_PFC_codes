function F_grav_B = gravity_body(m, R_IB, p)
% GRAVITY_BODY  Força gravitacional no referencial do corpo (B).
%
%   A gravidade no inercial é sempre g_I = [0; 0; -g] (Z_I aponta para cima).
%   Para usá-la na equação de Newton no referencial do corpo, deve-se
%   rotacionar para B usando a transposta de R_IB:
%
%     F_grav_B = m * R_IB' * g_I
%
%   Nota: R_IB' = R_BI transforma do inercial para o corpo.
%
% ENTRADAS
%   m    — massa atual [kg]
%   R_IB — matriz de rotação corpo→inercial (3×3)
%   p    — struct de parâmetros (contém p.g)
%
% SAÍDA
%   F_grav_B — força gravitacional no referencial do corpo [N]  (3×1)
%
% REFERÊNCIA: Tewari (2007), Aula 21 de MecVoo (Eq. de Newton no corpo)

% Vetor gravidade no inercial (Z_I aponta para cima → g em -Z_I)
g_I = [0; 0; -p.g];

% Rotacionar para o referencial do corpo
% R_IB' = R_BI: transforma de I para B
F_grav_B = m * (R_IB' * g_I);

end
