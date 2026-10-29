function [PC, L] = design_pid_budget_curve(gammas, nstart, maxeval, osTarget, tag)
%   osTarget (default 10): overshoot target of the SQP search; the final check
%   uses the specification 10 %% (as pid_estimator_yaw_margin; the SQP meets the active
%   OS constraint only to the solver tolerance). tag: output-file suffix.
%DESIGN_PID_BUDGET_CURVE PID candidates of the budget-performance curve
%   that fail the response requirement (gamma = 0.50: phi/theta settling
%   2.30 s; gamma = 0.75: yaw settling 68.8 s, OS 10.7 %), searched again with
%   the procedure of design_pid_strict (same subfunctions, copied from that file):
%   z and psi multistart SQP, phi/theta + xy joint multistart SQP, nstart
%   starts each, all constraints of the covariance design problem with the local increment ceilings
%   multiplied by gamma, the main cascade ceiling, and the response
%   requirement; starts around the stored curve candidate and the PID-Budget
%   class gains. (i) A general repair step (pid_repair,
%   the feasibility_repair procedure) is applied to every infeasible SQP result;
%   (ii) the gammas are processed in increasing order and the feasible design
%   of the previous gamma (feasible for the larger gamma, whose local ceilings
%   are looser) is included as a start and as a candidate, as in the original
%   budget-curve search ("each search includes its incumbent"). Writes
%   ../results/analysis/pid_curve.mat and pid_curve_verify.csv.
if nargin < 1, gammas = [0.5 0.75]; end
if nargin < 2, nstart = 20; end
if nargin < 3, maxeval = 900; end
if nargin < 4, osTarget = 10; end
if nargin < 5, tag = ''; end
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; Sspec = S; S.osmax = osTarget; cap = D.cap; req = [15 5 2 5]; rng(70004,'twister');
PC = struct('Gamma',{},'Controller',{}); L = table(); V = table(); ix = {[], 13:18, [], 31:36};
gammas = sort(gammas); inc = cell(1,4);   % feasible gains of the previous (smaller) gamma: feasible for a larger gamma too
lg = @(p) log(min(max(p,1e-4),1-1e-4)/(1-min(max(p,1e-4),1-1e-4)));
for gm = gammas
    j = find(strcmp({D.Curve.Family},'PID') & abs([D.Curve.Gamma] - gm) < 1e-9, 1); A0 = D.Curve(j).Controller; dub = gm*D.dub;
    co = A0; oa = [co.tauhat co.ndhat]; ok = true;
    for c = [2 4]
        g0 = A0.G(ix{c})'; starts = [{zeros(1,5)}, arrayfun(@(k) 0.6*randn(1,5), 2:nstart, 'uni', 0)]; res = cell(1,nstart);
        % two starts at the PID-Budget class gains (unscaled and scaled 0.9): the stored
        % curve candidate can be degenerate (yaw of gamma = 0.75)
        gp = D.C(1).G(ix{c})';
        starts{2} = [log(gp(1:4)./g0(1:4)), lg(gp(5)) - lg(g0(5))];
        starts{3} = [log(gp(1:4).*[.9^2 .9^3 .9 .9]./g0(1:4)), lg(gp(5)) - lg(g0(5))];
        if ~isempty(inc{c}), starts{4} = [log(inc{c}(1:4)./g0(1:4)), lg(inc{c}(5)) - lg(g0(5))]; end   % nested incumbent (smaller gamma)
        parfor k = 1:nstart, res{k} = fit_class(g0, starts{k}, c, S, dub, req, oa, maxeval); end
        for k = 1:nstart   % general repair step (pid_repair) for every infeasible SQP result
            if ~res{k}.feasible, [gr, okr, qr] = pid_repair(res{k}.g, c, S, [], dub(c), req, oa); if okr, res{k}.g = gr; res{k}.J = qr(1); res{k}.feasible = true; res{k}.repaired = true; end, end
        end
        if ~isempty(inc{c}), qi = channel_stats(2, inc{c}, c, S, [], oa, S.act, false, 4000); res{end+1} = struct('g', inc{c}, 'J', qi(1), 'feasible', true, 'exitflag', NaN, 'evals', 0); end
        [best, bj, nf] = pickbest(res, 'g'); L = [L; table(gm, c, nf, bj, 'VariableNames', {'Gamma','Class','Feasible','Objective'})]; %#ok<AGROW>
        if isempty(best), ok = false; break, end
        co.G(ix{c}) = best';
    end
    if ok
        b0 = [A0.G(19:24)', A0.G(1:6)']; starts = [{zeros(1,10)}, arrayfun(@(k) 0.6*randn(1,10), 2:nstart, 'uni', 0)]; res = cell(1,nstart);
        if ~isempty(inc{1}), starts{2} = [log(inc{3}(1:4)./b0(1:4)), lg(inc{3}(5)) - lg(b0(5)), log(inc{1}(1:4)./b0(7:10)), lg(inc{1}(5)) - lg(b0(11))]; end
        parfor k = 1:nstart, res{k} = fit_joint(b0, starts{k}, co, S, dub, req, cap, maxeval); end
        for k = 1:nstart   % general repair: inner class, then outer class with the repaired inner gains, then cascade check
            if res{k}.feasible, continue, end
            Gk = res{k}.G; [gi, oki] = pid_repair(Gk(19:24)', 3, S, [], dub(3), req, oa); if ~oki, continue, end
            [go, oko] = pid_repair(Gk(1:6)', 1, S, gi, dub(1), req, oa); if ~oko, continue, end
            cc = co; cc.G(19:24) = gi'; cc.G(25:30) = gi'; cc.G(1:6) = go'; cc.G(7:12) = go'; zc = cascade_stats(cc, S);
            if zc(4) <= cap*0.99, res{k}.G = cc.G; res{k}.J = zc(1); res{k}.feasible = true; end
        end
        if ~isempty(inc{1}), cc = co; cc.G(19:24) = inc{3}'; cc.G(25:30) = inc{3}'; cc.G(1:6) = inc{1}'; cc.G(7:12) = inc{1}'; zc = cascade_stats(cc, S); res{end+1} = struct('G', cc.G, 'J', zc(1), 'feasible', true, 'exitflag', NaN, 'evals', 0); end
        [best, bj, nf] = pickbest(res, 'G'); L = [L; table(gm, 13, nf, bj, 'VariableNames', {'Gamma','Class','Feasible','Objective'})]; %#ok<AGROW>
        if isempty(best), ok = false; else, co.G = best; end
    end
    if ok   % OS met only to the solver tolerance: reduce the set-point weight b (0.1 % steps), as in finalize_estimator_designs
        for h = [1 2 3 4]
            ii = {1:6, 13:18, 19:24, 31:36}; g = co.G(ii{h})'; gin = []; if h == 1, gin = co.G(19:24)'; end
            [os, ~] = step_response_metrics(chan_model(2, g, S.Ts, gin, S.act, [], oa), S.Ts); n = 0;
            while os > S.osmax && n < 200, g(5) = g(5)*0.999; n = n + 1; [os, ~] = step_response_metrics(chan_model(2, g, S.Ts, gin, S.act, [], oa), S.Ts); end
            co.G(ii{h}) = g'; if h == 1, co.G(7:12) = g'; elseif h == 3, co.G(25:30) = g'; end
        end
        co.name = sprintf('PID curve gamma=%.2f', gm); V = [V; verify_design(co, Sspec, dub, req, cap)];
        inc = {co.G(1:6)', co.G(13:18)', co.G(19:24)', co.G(31:36)'};
    else, co.name = sprintf('PID curve gamma=%.2f INFEASIBLE', gm); end %#ok<AGROW>
    PC(end+1) = struct('Gamma', gm, 'Controller', co); %#ok<AGROW>
end
save(fullfile(RA,['pid_curve' tag '.mat']), 'PC', 'L', 'V'); writetable(V, fullfile(RA,['pid_curve_verify' tag '.csv'])); writetable(L, fullfile(RA,['pid_curve_search' tag '.csv']));
end

function [best, bj, nf] = pickbest(res, fld)
best = []; bj = inf; nf = 0;
for k = 1:numel(res), r = res{k}; if r.feasible, nf = nf + 1; if r.J < bj, bj = r.J; best = r.(fld); end, end, end
end

% =====================================================================
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

function r = fit_joint(b0, x0, PS, S, dub, req, cap, maxeval)
oa = [PS.tauhat PS.ndhat];
    function [gi, go] = dec(x)
        gi = [b0(1:4).*exp(x(1:4)), 1/(1+exp(-(logit(b0(5)) + x(5)))), 0];
        go = [b0(7:10).*exp(x(6:9)), 1/(1+exp(-(logit(b0(11)) + x(10)))), 0];
    end
cache = containers.Map('KeyType','char','ValueType','any');
    function z = ev(x, nf)
        key = sprintf('%.12g,', [x nf]);
        if isKey(cache, key), z = cache(key); return, end
        [gi, go] = dec(x); co = PS; co.G(19:24) = gi'; co.G(25:30) = gi'; co.G(1:6) = go'; co.G(7:12) = go';
        qi = channel_stats(2, gi, 3, S, [], oa, S.act, false, nf);
        [osi, tsi] = step_response_metrics(chan_model(2, gi, S.Ts, [], S.act, [], oa), S.Ts);
        qo = channel_stats(2, go, 1, S, gi, oa, S.act, false, nf);
        [oso, tso] = step_response_metrics(chan_model(2, go, S.Ts, gi, S.act, [], oa), S.Ts);
        zc = cascade_stats(co, S);
        z = [zc(1), qi(2), qi(3), qi(4), osi, tsi, qi(6), qo(2), qo(3), qo(4), oso, tso, qo(6), zc(4)]; cache(key) = z;
    end
    function f = obj(x), z = ev(x, 600); f = z(1); if ~isfinite(f), f = 1e3; end, end
    function [v, e] = con(x)
        z = ev(x, 600);
        v = [z(2)/S.ub(3)-1, z(3)/dub(3)-1, z(4)/S.Msmax-1, z(5)/S.osmax-1, z(6)/req(3)-1, z(7)-0.9999, ...
             z(8)/S.ub(1)-1, z(9)/dub(1)-1, z(10)/S.Msmax-1, z(11)/S.osmax-1, z(12)/req(1)-1, z(13)-0.9999, z(14)/(cap*0.99)-1];
        v(~isfinite(v)) = 1e3; e = [];
    end
opts = optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',maxeval, ...
    'MaxIterations',80,'StepTolerance',1e-6,'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-7,'FiniteDifferenceStepSize',1e-4);
lb = -3*ones(1,10); ub = 3*ones(1,10); lb([5 10]) = -10; ub([5 10]) = 10;
try
    [x, ~, ef, out] = fmincon(@obj, max(lb, min(ub, x0)), [], [], [], [], lb, ub, @con, opts);
catch
    x = x0; ef = -99; out.funcCount = 0;
end
[gi, go] = dec(x); G = PS.G; G(19:24) = gi'; G(25:30) = gi'; G(1:6) = go'; G(7:12) = go';
r = struct('G', G, 'J', inf, 'feasible', false, 'exitflag', ef, 'evals', out.funcCount);
z = ev(x, 4000);
v = [z(2)/S.ub(3)-1, z(3)/dub(3)-1, z(4)/S.Msmax-1, z(5)/S.osmax-1, z(6)/req(3)-1, z(7)-0.9999, ...
     z(8)/S.ub(1)-1, z(9)/dub(1)-1, z(10)/S.Msmax-1, z(11)/S.osmax-1, z(12)/req(1)-1, z(13)-0.9999, z(14)/(cap*0.99)-1];
r.feasible = all(isfinite(z)) && all(v <= 1e-7);
if r.feasible, r.J = z(1); end
end

function y = logit(p), p = min(max(p, 1e-4), 1-1e-4); y = log(p/(1-p)); end
