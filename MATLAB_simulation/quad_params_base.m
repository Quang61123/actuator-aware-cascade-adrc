function pv = quad_params_base(scale)
%QUAD_PARAMS_BASE Nominal quadrotor parameter vector (Luukkonen 2011; also used in
%the reference model of Vo Van & Ha Duy 2025, Table 2).
%   pv = [m g l k b IM Ixx Iyy Izz A wmax Tmax tau_max(3)]'
%   scale (optional) = [s_m s_I s_A] multiplies mass, inertia and drag
%   (used for the model-mismatch scenario).
if nargin < 1, scale = [1 1 1]; end
m    = 0.468*scale(1);      % mass [kg]
g    = 9.81;                % gravity [m/s^2]
l    = 0.225;               % arm length [m]
k    = 2.980e-6;            % thrust coefficient [N s^2]
b    = 1.140e-7;            % drag-torque coefficient [N m s^2]
IM   = 3.357e-5;            % rotor inertia [kg m^2]
Ixx  = 4.856e-3*scale(2);   % [kg m^2]
Iyy  = 4.856e-3*scale(2);
Izz  = 8.801e-3*scale(2);
A    = 0.25*scale(3);       % linear drag [kg/s] (Luukkonen 2011)
wmax = 1000;                % rotor speed limit [rad/s] (assumed)
Tmax = 2*0.468*9.81;        % collective thrust limit [N] (assumed, T/W = 2)
taumax = [0.2 0.2 0.05];    % torque limits [N m] (assumed)
pv = [m g l k b IM Ixx Iyy Izz A wmax Tmax taumax]';
end


