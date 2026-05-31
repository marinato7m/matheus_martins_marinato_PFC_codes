function r_cg = cg_model(m, p)
% CG_MODEL  Posição do CG em relação à base do foguete (bocal).
%
%   r_cg = cg_model(m, p)
%
%   Modelo: interpolação linear entre os valores extremos.
%   O CG se desloca para BAIXO conforme o propelente queima:
%     - Motores Merlin (4.410 kg) fixos na base (~0,5 m)
%     - Propelente (tanques acima do CG) é consumido
%     → CG migra em direção à base
%
%   r_CG(m) = r_CG_f + (r_CG_0 - r_CG_f) * (m - mf) / (m0 - mf)
%
%   Esta função é chamada a cada passo de integração em rocket_6dof.m
%   para atualizar o braço de momento do TVC: L_TVC = r_cg.
%
%   As equações de Newton/Euler permanecem válidas pois são escritas
%   para o CG instantâneo (Tewari, 2007, seção 4.3).
%
% ENTRADA
%   m — massa atual [kg]
%   p — struct de parâmetros (falcon9_params)
%
% SAÍDA
%   r_cg — distância CG-bocal [m]
%          (= braço de alavanca do TVC = L_TVC)

frac  = (m - p.mf) / (p.m0 - p.mf);
frac  = max(0.0, min(1.0, frac));     % saturação: 0 ≤ frac ≤ 1
r_cg  = p.r_CG_f + (p.r_CG_0 - p.r_CG_f) * frac;

end
