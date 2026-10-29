function xdot = quad_dynamics_base(x, w, Fd, Td, pv)
%QUAD_DYNAMICS_BASE 6-DOF Newton-Euler quadrotor model (12 states).
%   x = [x y z vx vy vz phi theta psi p q r]'  (inertial pos/vel, ZYX Euler, body rates)
%   w = rotor speeds [rad/s]; Fd = inertial disturbance force [N];
%   Td = body disturbance torque [N m]; pv from quad_params_base.
%#codegen
m=pv(1); g=pv(2); l=pv(3); k=pv(4); b=pv(5); IM=pv(6);
I = pv(7:9); A = pv(10);
v = x(4:6); phi=x(7); th=x(8); psi=x(9); nu = x(10:12);
w2 = w(:).^2;
T   = k*sum(w2);
tau = [l*k*(w2(4)-w2(2)); l*k*(w2(3)-w2(1)); b*(-w2(1)+w2(2)-w2(3)+w2(4))];
cf=cos(phi); sf=sin(phi); ct=cos(th); st=sin(th); cp=cos(psi); sp=sin(psi);
R = [cp*ct, cp*st*sf-sp*cf, cp*st*cf+sp*sf;
     sp*ct, sp*st*sf+cp*cf, sp*st*cf-cp*sf;
     -st,   ct*sf,          ct*cf];
acc = [0;0;-g] + R*[0;0;T]/m - A*v/m + Fd(:)/m;
Om = w(1)-w(2)+w(3)-w(4);
p=nu(1); q=nu(2); r=nu(3);
nud = [((I(2)-I(3))*q*r - IM*q*Om)/I(1);
       ((I(3)-I(1))*p*r + IM*p*Om)/I(2);
        (I(1)-I(2))*p*q/I(3)] + (tau + Td(:))./I(:);
Wi = [1, sf*st/ct, cf*st/ct; 0, cf, -sf; 0, sf/ct, cf/ct];
xdot = [v; acc; Wi*nu; nud];
end


