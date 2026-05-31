function [F_tvc_B, M_tvc_B] = tvc_model(T, delta_y, delta_z, L_TVC)
% TVC_MODEL  Força e momento do TVC no referencial do corpo (B).
%
%   [F_tvc_B, M_tvc_B] = tvc_model(T, delta_y, delta_z, L_TVC)
%
% ── CONVENÇÃO DE SINAIS ──────────────────────────────────────────────────
%
%   delta_y > 0  →  nariz inclina para +Y_B
%   delta_z > 0  →  nariz inclina para +X_B
%
% ── DERIVAÇÃO DO SINAL CORRETO ───────────────────────────────────────────
%
%   O bocal está em r_bocal = [0; 0; -L_TVC] (abaixo do CG, pois Z_B
%   aponta para o nariz).
%
%   Para nariz ir para +Y_B:
%     - É necessário Mx < 0 (rotação ao redor de -X_B)
%     - Mx = cross(r_bocal, F_tvc)[1] = L_TVC * Fy
%     - Para Mx < 0 → Fy < 0
%     - Logo: Fy = -T * sin(delta_y)  ← sinal negativo é correto
%
%   Para nariz ir para +X_B:
%     - É necessário My < 0 (rotação ao redor de -Y_B)
%     - My = cross(r_bocal, F_tvc)[2] = -L_TVC * Fx
%     - Para My < 0 → Fx > 0
%     - Logo: Fx = +T * sin(delta_z)
%
% ── FORÇA DO TVC (pequenos ângulos) ──────────────────────────────────────
%
%   F_x =  T * sin(delta_z)    → gera My → controla theta (pitch)
%   F_y = -T * sin(delta_y)    → gera Mx → inclina para Y_B
%   F_z =  T * cos(delta_y) * cos(delta_z)   → componente principal
%
% ── MOMENTOS RESULTANTES ─────────────────────────────────────────────────
%
%   M = cross([0; 0; -L_TVC], F_tvc)
%     = [L_TVC * Fy;  -L_TVC * Fx;  0]
%
%   Mx = -L_TVC * T * sin(delta_y)  → <0 para delta_y>0 → nariz +Y_B ✓
%   My = -L_TVC * T * sin(delta_z)  → <0 para delta_z>0 → nariz +X_B ✓
%   Mz = 0                           → motor central sem controle de roll
%
% ── VERIFICAÇÃO RÁPIDA ───────────────────────────────────────────────────
%
%   [F, M] = tvc_model(845000, 0, 0, 18)
%     F deve ser [0; 0; 845000]  (empuxo puro em Z_B, sem deflexão)
%     M deve ser [0; 0; 0]       (sem momento)
%
%   [F, M] = tvc_model(845000, 0.1, 0, 18)
%     Mx deve ser ≈ -845000*18*0.1 = -1.521e6 N·m  (negativo → nariz +Y_B)
%
%   [F, M] = tvc_model(845000, 0, 0.1, 18)
%     My deve ser ≈ -845000*18*0.1 = -1.521e6 N·m  (negativo → nariz +X_B)
%
% ENTRADAS
%   T       — magnitude do empuxo [N]  (sempre positivo)
%   delta_y — deflexão TVC para Y [rad]  (>0 → nariz para +Y_B)
%   delta_z — deflexão TVC para X [rad]  (>0 → nariz para +X_B)
%   L_TVC   — braço de momento CG→bocal [m]  (output de cg_model)
%
% SAÍDAS
%   F_tvc_B — força do TVC no corpo [N]   (3×1)
%   M_tvc_B — momento do TVC no corpo [N·m] (3×1)
%
% REFERÊNCIA: Tewari (2007), cap. 8; Jenie et al. (2019)

% Força do TVC (pequenos ângulos — erro < 0.6% para |δ| < 6°)
F_tvc_B = [ T * sin(delta_z);       %  Fx → gera My → nariz +X_B
            -T * sin(delta_y);       %  Fy → gera Mx → nariz +Y_B (sinal negativo!)
             T * cos(delta_y) * cos(delta_z)];  % componente principal

% Posição do bocal em relação ao CG (Z_B aponta para cima → bocal ABAIXO)
r_bocal = [0; 0; -L_TVC];

% Momento: produto vetorial
% cross([0;0;-L], [Fx;Fy;Fz]) = [L*Fy; -L*Fx; 0]
M_tvc_B = cross(r_bocal, F_tvc_B);

% Resultado analítico (para conferência no debug):
%   Mx = -L_TVC * T * sin(delta_y)   → nariz para +Y_B quando Mx < 0 ✓
%   My = -L_TVC * T * sin(delta_z)   → nariz para +X_B quando My < 0 ✓
%   Mz = 0

end