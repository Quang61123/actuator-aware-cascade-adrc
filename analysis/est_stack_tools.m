function out = est_stack_tools(what, varargin)
%EST_STACK_TOOLS Estimator-stack evaluation used by the -E designs.
%   Delay convention, fixed before the design: the
%   filtered body rates (roll, pitch AND yaw) are received 6 ms late; position,
%   velocity and attitude estimates have no measurable lag. It is applied to
%   the sampled models used for Ms, the step response and the covariances.
%   Noise model: pooled estimator-coloured models of the estimator campaign (est_noise_models).
%     ds = est_stack_tools('ds', NF)                      delay vector (12 signals)
%     dl = est_stack_tools('dl', ds, c)                   [d_y d_v d_th d_om] of class c
%     M  = est_stack_tools('model', kind, g, c, gin, oa, ds)
%     ms = est_stack_tools('ms', kind, g, c, gin, oa, ds, nfreq)
%     [os ts] via est_stack_tools('resp', ...) -> out = [os ts]
%     z  = est_stack_tools('cascade', co, S, NF, ds)      [sigma_e sigma_u sigma_du sigma_du_inner rho]
%     z  = est_stack_tools('single', kind, g, c, S, NF, ds, oa)  [sigma_e sigma_u sigma_du ~ rho] (z, psi)
switch what
    case 'ds'
        NF = varargin{1}; out = NF.dsig; out(10:12) = 6;
    case 'dl'
        ds = varargin{1}; c = varargin{2};
        switch c
            case 1, out = ds([1 4 8 11]);
            case 2, out = [ds([3 6]) 0 0];
            case 3, out = [ds([8 11]) 0 0];
            case 4, out = [ds([9 12]) 0 0];
        end
    case 'model'
        [kind, g, c, gin, oa, ds] = varargin{:}; S = sset();
        out = chan_model_est(kind, g, S.Ts, gin, S.act, [], oa, est_stack_tools('dl', ds, c));
    case 'ms'
        [kind, g, c, gin, oa, ds, nf] = varargin{:}; S = sset();
        out = msens(est_stack_tools('model', kind, g, c, gin, oa, ds), S.Ts, nf);
    case 'resp'
        [kind, g, c, gin, oa, ds] = varargin{:}; S = sset();
        [os, ts] = step_response_metrics(est_stack_tools('model', kind, g, c, gin, oa, ds), S.Ts); out = [os ts];
    case 'cascade'
        [co, S, NF, ds] = varargin{:}; oa = [co.tauhat co.ndhat];
        MX = chan_model_est(co.kind, co.G(1:6)', S.Ts, co.G(19:24)', S.act, [], oa, ds([1 4 8 11]));
        MY = chan_model_est(co.kind, co.G(1:6)', S.Ts, co.G(19:24)', S.act, [], oa, ds([2 5 7 10]));
        zX = cov_estimator_aware(MX, NF.X_pooled); zY = cov_estimator_aware(MY, NF.Y_pooled);
        out = inf(1,5); out(5) = max(zX(5), zY(5)); if out(5) >= 1, return, end
        vd = rwpart(MX, S.qd(1), S.Ts); if ~isfinite(vd), return, end
        out = [sqrt(zX(1)^2 + vd), max(zX(2),zY(2)), max(zX(3),zY(3)), max(zX(4),zY(4)), out(5)];
    case 'single'
        [kind, g, c, S, NF, ds, oa] = varargin{:};
        M = chan_model_est(kind, g, S.Ts, [], S.act, [], oa, est_stack_tools('dl', ds, c));
        if c == 2, nm = NF.Z_pooled; else, nm = NF.P_pooled; end
        z = cov_estimator_aware(M, nm); out = inf(1,5); out(5) = z(5); if z(5) >= 1, return, end
        vd = rwpart(M, S.qd(c), S.Ts); if ~isfinite(vd), return, end
        out = [sqrt(z(1)^2 + vd), z(2), z(3), NaN, z(5)];
end
end

function S = sset()
persistent SS
if isempty(SS)
    here = fileparts(mfilename('fullpath'));
    D = load(fullfile(here,'..','results','design_gains.mat'), 'S'); SS = D.S;
end
S = SS;
end

function vd = rwpart(M, qd, Ts)
k = M.keep; F = M.A(k,k); bd = M.A(k,M.id); ce = M.Cy(k); xss = (eye(numel(k)) - F)\bd;
if abs(ce*xss) > 1e-7*max(1, norm(ce)*norm(xss)), vd = inf; return, end
vd = qd*Ts*(xss'*dlyap(F', ce'*ce)*xss);
end

function Ms = msens(M, Ts, nf)
% input-sensitivity peak at the virtual-input summing point (as chan_ms, refined)
k = M.keep; A = M.A(k,k); B = M.B(k,7); C = M.Ct(k); Dd = M.Dt(7);
if max(abs(eig(A))) >= 1, Ms = Inf; return, end
sys = ss(A,B,C,Dd,Ts); wg = [0 logspace(-2,log10(pi/Ts),nf)];
mag = abs(squeeze(freqresp(sys,wg))).'; Ms = max(mag);
pk = find(mag(2:end-1) >= mag(1:end-2) & mag(2:end-1) >= mag(3:end)) + 1; pk = pk(mag(pk) > 1 + 1e-8);
opts = optimset('TolX',1e-7,'Display','off');
for j = pk
    [~, neg] = fminbnd(@(lw) -abs(evalfr(sys, exp(1i*exp(lw)*Ts))), log(max(wg(j-1),1e-9)), log(wg(j+1)), opts);
    Ms = max(Ms, -neg);
end
end
