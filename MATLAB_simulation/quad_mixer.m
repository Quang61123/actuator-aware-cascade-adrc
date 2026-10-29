function [w, Tact, tauact] = quad_mixer(T, tau, pv)
%QUAD_MIXER Inverts the '+' rotor mixing (Luukkonen 2011) with saturation.
%#codegen
k = pv(4); l = pv(3); b = pv(5); wmax = pv(11);
T   = min(max(T,0),pv(12));
tau = min(max(tau(:),-pv(13:15)),pv(13:15));
M = [ k    k    k    k;
      0   -l*k  0    l*k;
     -l*k  0    l*k  0;
     -b    b   -b    b];
w2 = M\[T; tau];
w2 = min(max(w2,0),wmax^2);
w  = sqrt(w2);
F  = M*w2;
Tact = F(1); tauact = F(2:4);
end


