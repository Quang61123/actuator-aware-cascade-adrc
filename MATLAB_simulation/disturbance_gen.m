function [Fd, Td] = disturbance_gen(t, scen)
%DISTURBANCE_GEN Deterministic force/torque steps (S3, S4) and a torque
%   disturbance in the turbulence scenarios (S5, S6). Wind turbulence is an
%   exogenous signal (make_exo) added to Fd by the simulators.
%#codegen
Fd = zeros(3,1); Td = zeros(3,1);
if scen == 3 || scen == 4
    if t >= 10, Fd = [0.3; -0.3; -0.3]; end
    if t >= 20, Td = [0.01; -0.01; 0.005]; end
elseif scen == 5 || scen == 6
    if t >= 15, Td = [0.005*sin(2*pi*0.5*t); 0.005; 0]; end
end
end


