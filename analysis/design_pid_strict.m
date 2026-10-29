function [PS, L] = design_pid_strict(nstart, maxeval)
%DESIGN_PID_STRICT Strict PID design (PID-Budget-S) under the response requirement.
%   The stored PID-Budget-S fails the response requirement (yaw: 2% settling
%   68.8 s, overshoot 10.7 % on a 400-s horizon). It is searched again with
%   the same constraints as all other designs plus the response requirement:
%     every class: sigma_u <= ub(c), local sigma_du <= dub(c), Ms <= 1.98,
%       OS <= 10 % and 2% settling <= req(c) (unit step of the sampled linear
%       closed loop, 400-s horizon; xy on the cascade model), Schur stable;
%     strict cascade ceiling: inner roll increment of the full cascade
%       <= 0.598 (1% margin, as in the other strict designs).
%   Search (nstart >= 20 starts per search):
%     z and psi: fmincon (SQP) over the five PID gains (log scale, set-point
%       weight in logit form), objective sigma_e of the class, starts =
%       PID-Budget gains, their 0.9 scaling and random log offsets (sd 0.6);
%     phi/theta + xy jointly (cascade): fmincon over the ten inner/outer
%       gains, objective cascade sigma_e of xy, starts = stored
%       PID-Budget-S, PID-Budget (inner scaled 0.9) and random log offsets.
%   Every candidate is verified with Ms on 4000 frequencies plus refinement
%   and the exact 400-s step; the feasible candidate with the smallest
%   objective is kept. If none is feasible, PS is empty and this is reported.
%   Writes ../results/analysis/pidS.mat and pidS_starts.csv.
if nargin < 1, nstart = 20; end
if nargin < 2, maxeval = 900; end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub; cap = D.strict;
req = [15 5 2 5]; rng(70001, 'twister'); L = table();
P0 = D.C(1); PS = D.Strict(1); oa = [P0.tauhat P0.ndhat];
% ---------- z and psi ----------
ix = {[], 13:18, [], 31:36};
for c = [2 4]
    g0 = P0.G(ix{c})'; gs = g0; gs(1:4) = gs(1:4).*[.9^2 .9^3 .9 .9];
    starts = cell(1, nstart); starts{1} = zeros(1,5); starts{2} = [log(gs(1:4)./g0(1:4)) 0];
    for k = 3:nstart, starts{k} = [0.6*randn(1,4) 0.6*randn]; end
    res = cell(1, nstart);
    parfor k = 1:nstart
        res{k} = fit_class(g0, starts{k}, c, S, dub, req, oa, maxeval);
    end
    best = []; bj = inf;
    for k = 1:nstart
        r = res{k}; L = [L; table(c, k, r.J, r.feasible, r.exitflag, r.evals, 'VariableNames', {'Class','Start','Objective','Feasible','ExitFlag','Evaluations'})]; %#ok<AGROW>
        if r.feasible && r.J < bj, bj = r.J; best = r.g; end
    end
    fprintf('PIDS class %d: best sigma_e %.5g (%d/%d feasible)\n', c, bj, sum(cellfun(@(r) r.feasible, res)), nstart);
    if isempty(best), PS = []; save(fullfile(RA,'pidS.mat'),'PS','L'); writetable(L, fullfile(RA,'pidS_starts.csv')); return, end
    PS.G(ix{c}) = best';
end
% ---------- phi/theta + xy jointly ----------
b0 = [D.Strict(1).G(19:24)', D.Strict(1).G(1:6)'];      % base: stored PID-Budget-S
g9 = P0.G(19:24)'; g9(1:4) = g9(1:4).*[.9^2 .9^3 .9 .9];
starts = cell(1, nstart); starts{1} = zeros(1,10);
starts{2} = [log(g9(1:4)./b0(1:4)), logit(g9(5)) - logit(b0(5)), zeros(1,5)];
starts{3} = [log(P0.G(19:22)'./b0(1:4)), logit(P0.G(23)) - logit(b0(5)), log(P0.G(1:4)'./b0(7:10)), logit(P0.G(5)) - logit(b0(11))];
for k = 4:nstart, starts{k} = 0.6*randn(1,10); end
res = cell(1, nstart);
parfor k = 1:nstart
    res{k} = fit_joint(b0, starts{k}, PS, S, dub, req, cap, maxeval);
end
best = []; bj = inf;
for k = 1:nstart
    r = res{k}; L = [L; table(13, k, r.J, r.feasible, r.exitflag, r.evals, 'VariableNames', {'Class','Start','Objective','Feasible','ExitFlag','Evaluations'})]; %#ok<AGROW>
    if r.feasible && r.J < bj, bj = r.J; best = r.G; end
end
fprintf('PIDS joint: best sigma_e %.5g (%d/%d feasible)\n', bj, sum(cellfun(@(r) r.feasible, res)), nstart);
writetable(L, fullfile(RA,'pidS_starts.csv'));
if isempty(best), PS = []; save(fullfile(RA,'pidS.mat'),'PS','L'); return, end
PS.G = best; PS.name = 'PID-Budget-S';
PS.note = 'strict PID searched with the response requirement';
V = verify_design(PS, S, dub, req, cap);
save(fullfile(RA,'pidS.mat'), 'PS', 'L', 'V');
writetable(V, fullfile(RA,'pidS_verify.csv'));
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
