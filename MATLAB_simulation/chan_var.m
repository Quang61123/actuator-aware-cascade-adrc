function [vy_d, vy_n, vu_n, vdu_n, rho] = chan_var(kind, gj, Ts, sy, sv, qd, gin, syi, svi, act, oact)
%CHAN_VAR Exact stationary variances of a channel (see chan_model):
%   vy_d  tracking variance due to the random-walk disturbance
%   vy_n  tracking variance due to measurement noise
%   vu_n  variance of the noise-induced virtual control
%   vdu_n variance of its sample-to-sample increment u(k)-u(k-1)
%         (high-frequency actuator activity)
if nargin < 11, oact = []; end
M = chan_model(kind, gj, Ts, gin, act, [], oact);
k = M.keep; A = M.A(k,k); bd = M.A(k,M.id);
Bn = M.B(k,2:5); cy = M.Cy(k); cu = M.Cu(k); dun = M.Du(2:5);
rho = max(abs(eig(A)));
if rho >= 1, vy_d = inf; vy_n = inf; vu_n = inf; vdu_n = inf; return, end
n = numel(k); xss = (eye(n) - A) \ bd;
if qd>0 && abs(cy*xss)>1e-7*max(1,norm(cy)*norm(xss))
    vy_d=inf; % nonzero DC gain makes random-walk output nonstationary
else
    vy_d = qd*Ts*(xss'*dlyap(A', cy'*cy)*xss);
end
W = diag([sy^2 sv^2 syi^2 svi^2]);
P = dlyap(A, Bn*W*Bn');
vy_n = cy*P*cy';
vu_n = cu*P*cu' + dun*W*dun';
% du_k = cu(A-I)x_{k-1} + (cu*Bn - dun) n_{k-1} + dun n_k
h = cu*(A - eye(n)); q = cu*Bn - dun;
vdu_n = h*P*h' + q*W*q' + dun*W*dun';
end



