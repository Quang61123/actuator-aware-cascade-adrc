function xdot = quad_dynamics(x, wcmd, Fd, Td, pv)
%QUAD_DYNAMICS 6-DOF model + first-order rotor dynamics.
%   x = [x y z vx vy vz phi theta psi p q r w1 w2 w3 w4]' (16 states)
%   wcmd = commanded rotor speeds (from the mixer), tau_m = pv(16).
%#codegen
w = max(x(13:16), 0);
xdot = [quad_dynamics_base(x(1:12), w, Fd, Td, pv); (wcmd(:) - x(13:16))/pv(16)];
end


