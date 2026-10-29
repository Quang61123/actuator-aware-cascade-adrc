function out = test_ctrl_frames(ctrlId, scen, Tend)
%TEST_CTRL_FRAMES Offline check of sitl_ctrl_step: the plant is fed
%   back through the NED/FRD conversion with exact (noise-free, undelayed)
%   states. Verifies the frame conversions before PX4 code generation.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
T = load(fullfile(here,'sitl_gain_table.mat'));
clear sitl_ctrl_step
pv = quad_params(); Ts = 1e-3; N = round(Tend/Ts);
x = zeros(16,1); x(13:16) = sqrt(pv(1)*pv(2)/(4*pv(4)));
X = zeros(N,16);
for i = 1:N
    [pos, vel, q, gyro] = state_to_px4(x);
    u = sitl_ctrl_step(pos, vel, q, gyro, true, ctrlId, scen, T.GT, T.KT, T.TH);
    wc = double(u(1:4))'*pv(11);
    X(i,:) = x';
    [Fd, Td] = disturbance_gen((i-1)*Ts, scen);
    h = Ts;
    k1 = quad_dynamics(x, wc, Fd, Td, pv); k2 = quad_dynamics(x+h/2*k1, wc, Fd, Td, pv);
    k3 = quad_dynamics(x+h/2*k2, wc, Fd, Td, pv); k4 = quad_dynamics(x+h*k3, wc, Fd, Td, pv);
    x = x + h/6*(k1+2*k2+2*k3+k4);
end
out.final = X(end,[1 2 3 9]);
out.X = X;
end
