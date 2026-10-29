function z = predict_estimator_stats(co, NF, ds, which, id)
%PREDICT_ESTIMATOR_STATS Estimator-aware hover predictions of a design.
%   which = 'pooled' uses NF.<class>_pooled; 'design' uses the per-design
%   model of controller id (NF.<class>(k) with ctrls == id).
%   ds: delay convention (est_stack_tools('ds')).
%   z = struct with fields (all standard deviations, hover, noise only):
%     se_x, se_y, se_z, se_psi : tracking errors (noise part)
%     su_phi, su_theta, su_z, su_psi : virtual controls
%     sdu_phi, sdu_theta, sdu_z, sdu_psi : 1-ms increments of the virtual controls
%     J : rotor-command increment [rad/s] (linearized mixer, uncorrelated classes)
%   Ordering of the logged virtual controls: R.ua = [u_phi u_theta u_psi u_z].
if nargin < 5, id = NaN; end
S = sset(); oa = [co.tauhat co.ndhat]; Ts = S.Ts; act = S.act;
mX = pick(NF, 'X', which, id); mY = pick(NF, 'Y', which, id); mZ = pick(NF, 'Z', which, id); mP = pick(NF, 'P', which, id);
MX = chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, ds([1 4 8 11]));
MY = chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, ds([2 5 7 10]));
MZ = chan_model_est(co.kind, co.G(13:18)', Ts, [], act, [], oa, [ds([3 6]) 0 0]);
MP = chan_model_est(co.kind, co.G(31:36)', Ts, [], act, [], oa, [ds([9 12]) 0 0]);
zX = cov_estimator_aware(MX, mX); zY = cov_estimator_aware(MY, mY);
zZ = cov_estimator_aware(MZ, mZ); zP = cov_estimator_aware(MP, mP);
% inner virtual control of the cascade: theta loop from class X, phi loop from class Y
suI = @(M, m) innerstd(M, m);
[suTh] = suI(MX, mX); [suPh] = suI(MY, mY);
z = struct('se_x', zX(1), 'se_y', zY(1), 'se_z', zZ(1), 'se_psi', zP(1), ...
    'su_phi', suPh, 'su_theta', suTh, 'su_z', zZ(2), 'su_psi', zP(2), ...
    'sdu_phi', zY(4), 'sdu_theta', zX(4), 'sdu_z', zZ(3), 'sdu_psi', zP(3));
pv = quad_params(); k = pv(4); l = pv(3); b = pv(5); wh = sqrt(pv(1)*pv(2)/(4*k));
Mx = [k k k k; 0 -l*k 0 l*k; -l*k 0 l*k 0; -b b -b b]; Jm = diag(1/(2*wh)*ones(4,1))/Mx;
I3 = pv(7:9); v = [(pv(1)*z.sdu_z)^2, (I3(1)*z.sdu_phi)^2, (I3(2)*z.sdu_theta)^2, (I3(3)*z.sdu_psi)^2];
z.J = sqrt(mean((Jm.^2)*v(:)));
end

function s = innerstd(M, m)
% standard deviation of the inner virtual control (Ci, Di) under the shaping filter
p = size(m.Cf,1); k = M.keep; F = M.A(k,k); Gn = M.B(k,2:1+p); nf = size(m.Af,1); n = numel(k);
At = [F, Gn*m.Cf; zeros(nf,n), m.Af]; Bt = [Gn*m.Df; m.Bf];
if max(abs(eig(At))) >= 1, s = inf; return, end
P = dlyap(At, Bt*Bt'); ci = [M.Ci(k), M.Di(2:1+p)*m.Cf]; di = M.Di(2:1+p)*m.Df;
s = sqrt(max(0, ci*P*ci' + di*di'));
end

function m = pick(NF, c, which, id)
if strcmp(which, 'pooled'), m = NF.([c '_pooled']); return, end
L = NF.(c); m = [];
for j = 1:numel(L), if any(L(j).ctrls == id) && numel(L(j).ctrls) == 1, m = L(j); end, end
assert(~isempty(m), 'predict_estimator_stats: no per-design model for id %d class %s', id, c);
end

function S = sset()
persistent SS
if isempty(SS)
    here = fileparts(mfilename('fullpath'));
    D = load(fullfile(here,'..','results','design_gains.mat'), 'S'); SS = D.S;
end
S = SS;
end
