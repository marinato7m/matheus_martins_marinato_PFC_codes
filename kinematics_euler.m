function euler_dot = kinematics_euler(euler, omega_B)
% KINEMATICS_EULER  Taxa de variação dos ângulos de Euler.
%
%   Relação cinemática entre as taxas angulares no corpo [p, q, r]
%   e as derivadas dos ângulos de Euler [phi_dot, theta_dot, psi_dot]:
%
%   [phi_dot  ]   [1   sin(phi)*tan(theta)   cos(phi)*tan(theta)] [p]
%   [theta_dot] = [0   cos(phi)             -sin(phi)           ] [q]
%   [psi_dot  ]   [0   sin(phi)/cos(theta)   cos(phi)/cos(theta)] [r]
%
%   SINGULARIDADE: a matriz L é singular em theta = ±90° (gimbal lock).
%   Para pouso vertical com desvios < 20°, esta condição nunca ocorre.
%
% ENTRADAS
%   euler   — [phi; theta; psi]  ângulos de Euler [rad]
%   omega_B — [p; q; r]          taxas angulares no corpo [rad/s]
%
% SAÍDA
%   euler_dot — [phi_dot; theta_dot; psi_dot]  derivadas [rad/s]
%
% REFERÊNCIA: Tewari (2007); equação padrão de mecânica de voo

phi   = euler(1);
theta = euler(2);
% psi = euler(3);   % não entra na matriz L

cp  = cos(phi);
sp  = sin(phi);
ct  = cos(theta);
tt  = tan(theta);   % = sin(theta)/cos(theta)

% Matriz cinemática de Euler (L)
L = [1,   sp*tt,   cp*tt ;
     0,   cp,     -sp    ;
     0,   sp/ct,   cp/ct ];

euler_dot = L * omega_B;

end
