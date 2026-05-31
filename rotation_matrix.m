function R_IB = rotation_matrix(phi, theta, psi)
% ROTATION_MATRIX  Matriz de rotação do corpo (B) para o inercial (I).
%
%   R_IB = rotation_matrix(phi, theta, psi)
%
%   Sequência ZYX (yaw → pitch → roll):
%     R_IB = Rz(psi) · Ry(theta) · Rx(phi)
%
%   Uso:
%     v_I = R_IB   * v_B   (corpo  → inercial)
%     v_B = R_IB'  * v_I   (inercial → corpo, transposta = inversa)
%
%   Convenção de eixos:
%     Z_B = nariz do foguete (eixo de simetria, aponta para cima)
%     Foguete vertical = (phi=0, theta=0, psi=0) → R_IB = eye(3)
%
%   TESTES DE SANIDADE (rodar após qualquer alteração):
%     assert(norm(rotation_matrix(0,0,0) - eye(3)) < 1e-12)   % identidade
%     R = rotation_matrix(0, pi/6, 0);
%     assert(norm(R*R' - eye(3)) < 1e-12)                     % ortogonalidade
%     nariz_I = R * [0;0;1]; assert(nariz_I(2) < 1e-12)       % pitch em YZ
%
% ENTRADAS
%   phi   — ângulo de roll  (rotação em torno de X_B) [rad]
%   theta — ângulo de pitch (rotação em torno de Y_B) [rad]
%   psi   — ângulo de yaw   (rotação em torno de Z_I) [rad]
%
% SAÍDA
%   R_IB  — matriz 3×3 ortogonal (det = +1)
%
% REFERÊNCIA: Tewari (2007), cap. 6; Jenie et al. (2019)

cphi = cos(phi);    sphi = sin(phi);
cth  = cos(theta);  sth  = sin(theta);
cpsi = cos(psi);    spsi = sin(psi);

%       [        col 1              col 2                    col 3        ]
R_IB = [cpsi*cth,   cpsi*sth*sphi - spsi*cphi,   cpsi*sth*cphi + spsi*sphi;
        spsi*cth,   spsi*sth*sphi + cpsi*cphi,   spsi*sth*cphi - cpsi*sphi;
        -sth,        cth*sphi,                     cth*cphi               ];

end
