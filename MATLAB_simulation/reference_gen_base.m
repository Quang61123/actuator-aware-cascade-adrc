function r = reference_gen_base(t, scen)
%REFERENCE_GEN_BASE r = [x y z psi xd yd zd xdd ydd zdd]'
%   scen < 5 : set-point (0.5, 0.5, 2) m, psi = 0.5 rad (step at t = 0)
%   scen = 5 : smooth climb to 2 m (0-5 s) then lemniscate (A = 1 m, period 20 s)
%#codegen
r = zeros(10,1);
if scen < 5
    r(1:4) = [0.5; 0.5; 2; 0.5];
else
    Tc = 5;
    if t < Tc
        s = t/Tc;
        r(3) = 2*(3*s^2 - 2*s^3); r(7) = 2*6*s*(1-s)/Tc; r(10) = 2*6*(1-2*s)/Tc^2;
    else
        r(3) = 2;
    end
    if t >= Tc
        tt = t - Tc; Om = 2*pi/20; Aa = 1;
        r(1) = Aa*sin(Om*tt);        r(5) = Aa*Om*cos(Om*tt);     r(8) = -Aa*Om^2*sin(Om*tt);
        r(2) = Aa/2*sin(2*Om*tt);    r(6) = Aa*Om*cos(2*Om*tt);   r(9) = -2*Aa*Om^2*sin(2*Om*tt);
    end
end
end


