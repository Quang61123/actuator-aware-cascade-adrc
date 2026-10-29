function L=observer_half_gain(wo,Ts)
% Predictive rate/disturbance ESO: exact mapped half-gain poles.
% Adapted from Herbst & Hempel, CEP 164 (2025) 106501, N=1.
% This is a rate-aided predictive adaptation, not their full current observer.
% L = Ad*Lc, Ad=[1 Ts;0 1], Lc=[1-z^2;(1-2*z*cos(h)+z^2)/Ts].
h=wo*Ts/2;z=exp(-h);
L=[2-2*z*cos(h);((1-z)^2+2*z*(1-cos(h)))/Ts];
end
