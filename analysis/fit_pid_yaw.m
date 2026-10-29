function r = fit_pid_yaw(g0, x0, S, dub, NF, ds, oa, Bpsi, maxeval)
%FIT_PID_YAW Standalone copy of the PID yaw search of design_estimator_stack (same code).
r = yaw_search(g0, x0, S, dub, NF, ds, oa, Bpsi, maxeval);
end

function r = yaw_search(g0, x0, S, dub, NF, ds, oa, Bpsi, maxeval)
dec = @(x) [g0(1:4).*exp(x(1:4)), 1/(1+exp(-(logit(g0(5)) + x(5)))), 0];
cache = containers.Map('KeyType','char','ValueType','any');
    function z = ev(x, nf)
        key = sprintf('%.12g,', [x nf]); if isKey(cache, key), z = cache(key); return, end
        g = dec(x); [~,~,vu,vdu] = chan_var(2, g, S.Ts, S.sy(4), S.sv(4), 0, [], 0, 0, S.act, oa);
        ze = est_stack_tools('single', 2, g, 4, S, NF, ds, oa); ms = est_stack_tools('ms', 2, g, 4, [], oa, ds, nf);
        rs = est_stack_tools('resp', 2, g, 4, [], oa, ds);
        z = [ze(1), sqrt(vu), sqrt(vdu), ze(3), ms, rs, ze(5)]; cache(key) = z;
    end
    function v = cv(z)
        v = [z(2)/S.ub(4)-1, z(3)/dub(4)-1, z(4)/(0.99*Bpsi)-1, z(5)/S.Msmax-1, z(6)/S.osmax-1, z(7)/5-1, z(8)-(1-1e-9)]; v(~isfinite(v)) = 1e3;   % Schur stability of the augmented (noise-filter) model
    end
    function f = obj(x), z = ev(x, 600); f = z(1); if ~isfinite(f), f = 1e3; end, end
    function [v, e] = con(x), v = cv(ev(x, 600)); e = []; end
opts = optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',maxeval, ...
    'MaxIterations',80,'StepTolerance',1e-6,'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-7,'FiniteDifferenceStepSize',1e-4);
lb = [-3 -3 -3 -3 -10]; ub = -lb;
try, [x, ~, ef, out] = fmincon(@obj, max(lb, min(ub, x0)), [], [], [], [], lb, ub, @con, opts); catch, x = x0; ef = -99; out.funcCount = 0; end
z = ev(x, 4000); r = struct('g', dec(x), 'J', inf, 'feasible', all(cv(z) <= 1e-7) && all(isfinite(z)), 'exitflag', ef, 'evals', out.funcCount);
if r.feasible, r.J = z(1); end
end

function y = logit(p), p = min(max(p, 1e-4), 1-1e-4); y = log(p/(1-p)); end
