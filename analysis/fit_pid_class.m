function r = fit_pid_class(g0, x0, c, S, dub, req, oa, maxeval)
%FIT_PID_CLASS Standalone copy of the per-class PID search of design_pid_strict (same code).
r = fit_class(g0, x0, c, S, dub, req, oa, maxeval);
end

function r = fit_class(g0, x0, c, S, dub, req, oa, maxeval)
dec = @(x) [g0(1:4).*exp(x(1:4)), 1/(1+exp(-(logit(g0(5)) + x(5)))), 0];
cache = containers.Map('KeyType','char','ValueType','any');
    function z = ev(x)
        key = sprintf('%.12g,', x);
        if isKey(cache, key), z = cache(key); return, end
        g = dec(x); q = channel_stats(2, g, c, S, [], oa, S.act, false, 600);
        [os, ts] = step_response_metrics(chan_model(2, g, S.Ts, [], S.act, [], oa), S.Ts);
        z = [q(1:4) os ts q(6)]; cache(key) = z;
    end
    function f = obj(x), z = ev(x); f = z(1); if ~isfinite(f), f = 1e3; end, end
    function [v, e] = con(x)
        z = ev(x); v = [z(2)/S.ub(c)-1, z(3)/dub(c)-1, z(4)/S.Msmax-1, z(5)/S.osmax-1, z(6)/req(c)-1, z(7)-0.9999];
        v(~isfinite(v)) = 1e3; e = [];
    end
opts = optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',maxeval, ...
    'MaxIterations',80,'StepTolerance',1e-6,'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-7,'FiniteDifferenceStepSize',1e-4);
lb = -3*ones(1,5); ub = 3*ones(1,5); lb(5) = -10; ub(5) = 10;
try
    [x, ~, ef, out] = fmincon(@obj, max(lb, min(ub, x0)), [], [], [], [], lb, ub, @con, opts);
catch
    x = x0; ef = -99; out.funcCount = 0;
end
g = dec(x); r = struct('g', g, 'J', inf, 'feasible', false, 'exitflag', ef, 'evals', out.funcCount);
q = channel_stats(2, g, c, S, [], oa, S.act, false, 4000);
[os, ts] = step_response_metrics(chan_model(2, g, S.Ts, [], S.act, [], oa), S.Ts);
r.feasible = all(isfinite(q(1:4))) && q(2) <= S.ub(c)*(1+1e-9) && q(3) <= dub(c)*(1+1e-9) && q(4) <= S.Msmax && q(6) < 1 && os <= S.osmax && ts <= req(c);
if r.feasible, r.J = q(1); end
end

function y = logit(p), p = min(max(p, 1e-4), 1-1e-4); y = log(p/(1-p)); end
