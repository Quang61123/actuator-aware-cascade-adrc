function [pos, vel, q, gyro, Rn] = state_to_px4(x)
%STATE_TO_PX4 Model state (x north, y west, z up; FLU body; ZYX Euler)
%   to PX4 conventions (NED world, FRD body). q = [w x y z] FRD->NED,
%   Rn = rotation matrix FRD->NED.
phi = x(7); th = -x(8); psi = -x(9);
pos = [x(1); -x(2); -x(3)]; vel = [x(4); -x(5); -x(6)];
gyro = [x(10); -x(11); -x(12)];
cr = cos(phi/2); sr = sin(phi/2); cp = cos(th/2); sp = sin(th/2); cy = cos(psi/2); sy = sin(psi/2);
q = [cr*cp*cy + sr*sp*sy; sr*cp*cy - cr*sp*sy; cr*sp*cy + sr*cp*sy; cr*cp*sy - sr*sp*cy];
if nargout > 4
    cf=cos(phi); sf=sin(phi); ct=cos(th); st=sin(th); cps=cos(psi); sps=sin(psi);
    Rn = [cps*ct, cps*st*sf-sps*cf, cps*st*cf+sps*sf;
          sps*ct, sps*st*sf+cps*cf, sps*st*cf-cps*sf;
          -st,    ct*sf,            ct*cf];
end
end
