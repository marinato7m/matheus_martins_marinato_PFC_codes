function F_aero_B = aero_model(vel_B, R_IB, p)
% AERO_MODEL  Força aerodinâmica no referencial do corpo (B).
%
%   F_aero_B = aero_model(vel_B, R_IB, p)
%
% ── HIPÓTESES ────────────────────────────────────────────────────────────
%
%   1. Apenas ARRASTO: D = 0.5 * rho * v_rel² * Sref * CD0
%      Sustentação e força lateral = 0 (ângulos de ataque pequenos, fase terminal)
%      Justificativa: Aula 24 do curso de mecânica de voo; Tewari (2007), cap. 7.
%
%   2. VENTO: p.v_vento = [0;0;0] no modo nominal.
%      Fácil de ativar — só mudar o parâmetro, sem alterar este arquivo.
%
%   3. DENSIDADE CONSTANTE: rho = 1.225 kg/m³ (baixa variação nas altitudes consideradas)
%
% ── FÍSICA ───────────────────────────────────────────────────────────────
%
%   Velocidade relativa ao vento (no corpo):
%     v_vento_B = R_IB' * v_vento_I
%     v_rel_B   = vel_B - v_vento_B
%
%   Pressão dinâmica:
%     q_din = 0.5 * rho * |v_rel_B|^2
%
%   Força de arrasto (opõe ao vetor de velocidade relativa):
%     F_aero_B = -D * v_rel_unit = -(q_din * Sref * CD0) * v_rel_B/|v_rel_B|
%
% ENTRADAS
%   vel_B — velocidade do CG no referencial do corpo [m/s]  (3×1)
%   R_IB  — matriz de rotação corpo→inercial (output de rotation_matrix)
%   p     — struct de parâmetros (falcon9_params)
%
% SAÍDA
%   F_aero_B — força aerodinâmica no referencial do corpo [N]  (3×1)
%
% REFERÊNCIA: Tewari (2007), cap. 7; Anderson (2015), cap. 5

% ── Velocidade relativa ao vento, no referencial do corpo ─────────────────
v_vento_I = p.v_vento;               % [0;0;0] nominal;
v_vento_B = R_IB' * v_vento_I;       % rotaciona para o corpo
v_rel_B   = vel_B - v_vento_B;        % velocidade aerodinâmica no corpo

% ── Magnitude ─────────────────────────────────────────────────────────────
v_mag = norm(v_rel_B);

% ── Guarda contra velocidade nula (evita divisão por zero) ────────────────
if v_mag < 1e-6
    F_aero_B = [0; 0; 0];
    return;
end

% ── Pressão dinâmica e arrasto ────────────────────────────────────────────
q_din    = 0.5 * p.rho * v_mag^2;    % [Pa]
D        = q_din * p.Sref * p.CD0;   % [N]

% Força opõe ao vetor de velocidade relativa
v_unit   = v_rel_B / v_mag;
F_aero_B = -D * v_unit;              % [N]  sempre aponta contra o movimento

end
