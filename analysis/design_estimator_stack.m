function DE = design_estimator_stack(nstart, maxeval)
%DESIGN_ESTIMATOR_STACK Designs for the PX4 + EKF2 stack.
%   Frozen references (written to frozen_ceilings.csv before the
%   search and not changed): B_cas^E = estimator-aware cascade roll/pitch
%   increment of PID-Budget and
%   B_psi^E = estimator-aware yaw increment of PID-Budget (value of the yaw-constrained
%   designs, 4.75e-4), both under the pooled noise model.
%   Delay convention: filtered body rates 6 ms, other signals no lag;
%   used for Ms, step response and covariances (est_stack_tools).
%   Specification for PID-E, AA-Budget-E and AA-HalfGain-E, all classes:
%     Ms <= 1.98 with the delay convention; OS <= 10 % and 2% settling
%     <= [15 5 2 5] s (400 s); local white-noise ceilings;
%     wo*Ts <= 0.2; cascade roll/pitch increment <= 0.99 B_cas^E and yaw
%     increment <= 0.99 B_psi^E (pooled noise model); no altitude ceiling.
%   A class of the primary design that already meets the specification is
%   kept (z for all three); the others are searched:
%     ADRC: yaw on a 20 x 40 grid of (wo/wc, wc); roll/pitch on a 12 x 12
%       grid of (wc, wo/wc), each feasible inner pair followed by the
%       horizontal search (20 ratios, 36-point scan plus bisection, as
%       cascade_outer_design_E); cost = estimator-aware tracking SD.
%     PID: multistart SQP (nstart starts) over the five gains (yaw) and the
%       ten roll/pitch + x/y gains jointly (cascade), same constraints.
%   Final check: Ms on 4000 frequencies plus refinement, exact step.
%   Writes ../results/analysis/designs_E.mat, designs_E_verify.csv,
%   design_E_search.csv, frozen_ceilings.csv.
if nargin < 1, nstart = 20; end
if nargin < 2, maxeval = 900; end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub;
L = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF = L.NF; ds = est_stack_tools('ds', NF);
E4 = load(fullfile(RP,'estimator_noise','est_designs_E.mat'), 'capE'); E5 = load(fullfile(RP,'yaw_ceiling','yaw_ceiling_designs.mat'), 'Bpsi');
Bcas = E4.capE; Bpsi = E5.Bpsi;
writetable(table(Bcas, Bpsi, {mat2str(ds)}, 'VariableNames', {'B_cas_E','B_psi_E','DelayConvention_samples'}), fullfile(RA,'frozen_ceilings.csv'));
req = [15 5 2 5]; rng(70002, 'twister'); search = table(); V = table();
base = D.C([1 2 4]); [base.note] = deal(''); names = {'PID-E','AA-Budget-E','AA-HalfGain-E'}; DE = base;
ck = fullfile(RA,'design_E_checkpoint.mat');
for f = 1:3
    if exist(ck,'file')
        K = load(ck); if numel(K.done) >= f && K.done(f), DE(f) = K.DE(f); V = K.V; search = K.search; continue, end
    end
    co = base(f); co.name = names{f}; oa = [co.tauhat co.ndhat];
    if co.kind == 2
        % ---- PID yaw (multistart) ----
        g0 = co.G(31:36)';
        if ~classok(co, 4, S, dub, NF, ds, Bcas, Bpsi)
            starts = [{zeros(1,5)}, arrayfun(@(k) 0.6*randn(1,5), 2:nstart, 'uni', 0)]; res = cell(1,nstart);
            parfor k = 1:nstart, res{k} = fit_pid_yaw(g0, starts{k}, S, dub, NF, ds, oa, Bpsi, maxeval); end
            [co, search] = pick(co, res, 31:36, search, names{f}, 'psi');
        end
        % ---- PID roll/pitch + x/y jointly ----
        b0 = [co.G(19:24)', co.G(1:6)'];
        starts = cell(1,nstart); starts{1} = zeros(1,10);
        starts{2} = [log([.9^2 .9^3 .9 .9]) 0 zeros(1,5)]; starts{3} = [log([.8^2 .8^3 .8 .8]) 0 log([.8^2 .8^3 .8 .8]) 0];
        for k = 4:nstart, starts{k} = 0.6*randn(1,10); end
        res = cell(1,nstart);
        parfor k = 1:nstart, res{k} = fit_pid_joint(b0, starts{k}, co, S, dub, NF, ds, Bcas, maxeval); end
        [co, search] = pick(co, res, [], search, names{f}, 'phth+xy');
    else
        % ---- ADRC yaw grid ----
        if ~classok(co, 4, S, dub, NF, ds, Bcas, Bpsi)
            rg = logspace(log10(1.5), log10(60), 20); rows = cell(1, numel(rg));
            parfor q = 1:numel(rg)
                r = rg(q); wg = logspace(log10(0.1), log10(S.wo_max/r), 40); R = nan(numel(wg), 6);
                for i = 1:numel(wg)
                    g = gains_of(co.kind, r, wg(i));
                    [ok, J, zz] = yawok(co.kind, g, S, dub, NF, ds, oa, Bpsi, 600);
                    R(i,:) = [r wg(i) J zz(1) zz(2) ok];
                end
                rows{q} = R;
            end
            R = vertcat(rows{:}); fe = R(R(:,6) == 1, :); [~, o] = sort(fe(:,3)); fe = fe(o,:); got = false;
            for i = 1:size(fe,1)                            % re-check best candidates on the fine grid
                g = gains_of(co.kind, fe(i,1), fe(i,2));
                if yawok(co.kind, g, S, dub, NF, ds, oa, Bpsi, 4000), co.G(31:36) = g'; got = true; break, end
            end
            search = [search; table(string(names{f}), "psi", size(R,1), size(fe,1), got, 'VariableNames', {'Design','Class','Evaluated','Feasible','Selected'})]; %#ok<AGROW>
            if ~got, co.name = [co.name '-INFEASIBLE']; DE(f) = co; continue, end
        end
        % ---- ADRC roll/pitch grid + horizontal search ----
        wcs = logspace(0, log10(co.G(19)), 12); rts = logspace(log10(1.5), log10(60), 12);
        [WC, RT] = ndgrid(wcs, rts); WC = WC(:); RT = RT(:); cand = cell(numel(WC),1); cost = inf(numel(WC),1);
        parfor q = 1:numel(WC)
            [cand{q}, cost(q)] = inner_outer(co, RT(q), WC(q), S, dub, NF, ds, Bcas);
        end
        [~, o] = sort(cost); got = false;
        for i = o(:)'
            if ~isfinite(cost(i)), break, end
            Vi = verify_design_stack(cand{i}, S, dub, NF, ds, Bcas, Bpsi, 4000);
            if all(Vi.Pass([1 3])), co = cand{i}; got = true; break, end
        end
        search = [search; table(string(names{f}), "phth+xy", numel(WC), sum(isfinite(cost)), got, 'VariableNames', {'Design','Class','Evaluated','Feasible','Selected'})]; %#ok<AGROW>
        if ~got, co.name = [co.name '-INFEASIBLE']; DE(f) = co; continue, end
    end
    co.note = 'designed for the PX4+EKF2 stack (delay convention, frozen ceilings)';
    DE(f) = co;
    V = [V; verify_design_stack(co, S, dub, NF, ds, Bcas, Bpsi, 4000)]; %#ok<AGROW>
    done = false(1,3); if exist(ck,'file'), K = load(ck); done(1:numel(K.done)) = K.done; end; done(f) = true;
    save(ck, 'DE', 'V', 'search', 'done');
    fprintf('EDESIGN %s done\n', co.name);
end
save(fullfile(RA,'designs_E.mat'), 'DE', 'V', 'search', 'Bcas', 'Bpsi', 'ds');
writetable(V, fullfile(RA,'designs_E_verify.csv')); writetable(search, fullfile(RA,'design_E_search.csv'));
end

function ok = classok(c0, h, S, dub, NF, ds, Bcas, Bpsi)
Vh = verify_design_stack(c0, S, dub, NF, ds, Bcas, Bpsi, 4000); ok = Vh.Pass(h);
end

% =====================================================================
function [co, search] = pick(co, res, ix, search, name, cls)
best = []; bj = inf; nf = 0;
for k = 1:numel(res)
    if res{k}.feasible, nf = nf + 1; if res{k}.J < bj, bj = res{k}.J; best = res{k}; end, end
end
search = [search; table(string(name), string(cls), numel(res), nf, ~isempty(best), 'VariableNames', {'Design','Class','Evaluated','Feasible','Selected'})];
if isempty(best), co.name = [co.name '-INFEASIBLE(' cls ')']; return, end
if isempty(ix), co.G = best.G; else, co.G(ix) = best.g'; end
end

function [ok, J, zz] = yawok(kind, g, S, dub, NF, ds, oa, Bpsi, nf)
J = inf; zz = [NaN NaN]; ok = false;
[~,~,vu,vdu] = chan_var(kind, g, S.Ts, S.sy(4), S.sv(4), 0, [], 0, 0, S.act, oa);
if ~(isfinite(vu) && sqrt(vu) <= S.ub(4)*(1+1e-9) && sqrt(vdu) <= dub(4)*(1+1e-9)), return, end
z = est_stack_tools('single', kind, g, 4, S, NF, ds, oa); zz = [z(3) NaN];
if ~(z(5) < 1 && z(3) <= 0.99*Bpsi), return, end
ms = est_stack_tools('ms', kind, g, 4, [], oa, ds, nf); zz(2) = ms; if ms > S.Msmax, return, end
rs = est_stack_tools('resp', kind, g, 4, [], oa, ds); if rs(1) > S.osmax || rs(2) > 5, return, end
ok = true; J = z(1);
end

function [best, bestJ] = inner_outer(co, rt, wc, S, dub, NF, ds, Bcas)
best = []; bestJ = inf; oa = [co.tauhat co.ndhat];
wo = rt*wc; if wo > S.wo_max, return, end
gi = gains_of(co.kind, rt, wc);
[~,~,vu,vdu] = chan_var(co.kind, gi, S.Ts, S.sy(3), S.sv(3), 0, [], 0, 0, S.act, oa);
if ~(isfinite(vu) && sqrt(vu) <= S.ub(3)*(1+1e-9) && sqrt(vdu) <= dub(3)*(1+1e-9)), return, end
if est_stack_tools('ms', co.kind, gi, 3, [], oa, ds, 600) > S.Msmax, return, end
rs = est_stack_tools('resp', co.kind, gi, 3, [], oa, ds); if rs(1) > S.osmax || rs(2) > 2, return, end
c = co; c.G(19:24) = gi'; c.G(25:30) = gi';
c0 = c; c0.G(1:6) = gains_of(c0.kind, 10, 0.05)'; c0.G(7:12) = c0.G(1:6);
z0 = est_stack_tools('cascade', c0, S, NF, ds); if ~(z0(4) <= 0.99*Bcas), return, end
ratios = logspace(log10(1.5), log10(60), 20);
for j = 1:numel(ratios)
    shape = ratios(j); grid = logspace(log10(.05), log10(S.wo_max/shape), 36); bw = NaN; idx = 0; cc = c;
    for i = numel(grid):-1:1
        cc.G(1:6) = gains_of(cc.kind, shape, grid(i))'; cc.G(7:12) = cc.G(1:6);
        if outer_ok(cc), bw = grid(i); idx = i; break, end
    end
    if ~isfinite(bw), continue, end
    if idx < numel(grid)
        lo = log(bw); hi = log(grid(idx+1));
        for b = 1:14
            w = exp((lo+hi)/2); cc.G(1:6) = gains_of(cc.kind, shape, w)'; cc.G(7:12) = cc.G(1:6);
            if outer_ok(cc), lo = log(w); else, hi = log(w); end
        end
        bw = exp(lo);
    end
    cc.G(1:6) = gains_of(cc.kind, shape, bw)'; cc.G(7:12) = cc.G(1:6);
    if outer_ok(cc)
        z = est_stack_tools('cascade', cc, S, NF, ds);
        if z(1) < bestJ, bestJ = z(1); best = cc; end
    end
end
    function ok = outer_ok(x)
        zw = cascade_stats(x, S);
        ok = all(isfinite(zw)) && zw(2) <= S.ub(1)*(1+1e-9) && zw(3) <= dub(1)*(1+1e-9); if ~ok, return, end
        ze = est_stack_tools('cascade', x, S, NF, ds); ok = all(isfinite(ze)) && ze(4) <= 0.99*Bcas; if ~ok, return, end
        ok = est_stack_tools('ms', x.kind, x.G(1:6)', 1, x.G(19:24)', oa, ds, 600) <= S.Msmax; if ~ok, return, end
        r = est_stack_tools('resp', x.kind, x.G(1:6)', 1, x.G(19:24)', oa, ds); ok = r(1) <= S.osmax && r(2) <= 15;
    end
end

function r = fit_pid_yaw(g0, x0, S, dub, NF, ds, oa, Bpsi, maxeval)
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

function r = fit_pid_joint(b0, x0, co, S, dub, NF, ds, Bcas, maxeval)
oa = [co.tauhat co.ndhat];
    function G = dec(x)
        gi = [b0(1:4).*exp(x(1:4)), 1/(1+exp(-(logit(b0(5)) + x(5)))), 0];
        go = [b0(7:10).*exp(x(6:9)), 1/(1+exp(-(logit(b0(11)) + x(10)))), 0];
        G = co.G; G(19:24) = gi'; G(25:30) = gi'; G(1:6) = go'; G(7:12) = go';
    end
cache = containers.Map('KeyType','char','ValueType','any');
    function z = ev(x, nf)
        key = sprintf('%.12g,', [x nf]); if isKey(cache, key), z = cache(key); return, end
        G = dec(x); c = co; c.G = G; gi = G(19:24)'; go = G(1:6)';
        [~,~,vui,vdui] = chan_var(2, gi, S.Ts, S.sy(3), S.sv(3), 0, [], 0, 0, S.act, oa);
        zw = cascade_stats(c, S); ze = est_stack_tools('cascade', c, S, NF, ds);
        msi = est_stack_tools('ms', 2, gi, 3, [], oa, ds, nf); ri = est_stack_tools('resp', 2, gi, 3, [], oa, ds);
        mso = est_stack_tools('ms', 2, go, 1, gi, oa, ds, nf); ro = est_stack_tools('resp', 2, go, 1, gi, oa, ds);
        z = [ze(1), sqrt(vui), sqrt(vdui), msi, ri, zw(2), zw(3), mso, ro, ze(4), max(ze(5), zw(5))]; cache(key) = z;
    end
    function v = cv(z)
        v = [z(2)/S.ub(3)-1, z(3)/dub(3)-1, z(4)/S.Msmax-1, z(5)/S.osmax-1, z(6)/2-1, ...
             z(7)/S.ub(1)-1, z(8)/dub(1)-1, z(9)/S.Msmax-1, z(10)/S.osmax-1, z(11)/15-1, z(12)/(0.99*Bcas)-1, z(13)-(1-1e-9)];
        v(~isfinite(v)) = 1e3;
    end
    function f = obj(x), z = ev(x, 600); f = z(1); if ~isfinite(f), f = 1e3; end, end
    function [v, e] = con(x), v = cv(ev(x, 600)); e = []; end
opts = optimoptions('fmincon','Algorithm','sqp','Display','off','MaxFunctionEvaluations',maxeval, ...
    'MaxIterations',80,'StepTolerance',1e-6,'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-7,'FiniteDifferenceStepSize',1e-4);
lb = -3*ones(1,10); ub = -lb; lb([5 10]) = -10; ub([5 10]) = 10;
try, [x, ~, ef, out] = fmincon(@obj, max(lb, min(ub, x0)), [], [], [], [], lb, ub, @con, opts); catch, x = x0; ef = -99; out.funcCount = 0; end
z = ev(x, 4000); r = struct('G', dec(x), 'J', inf, 'feasible', all(cv(z) <= 1e-7) && all(isfinite(z)), 'exitflag', ef, 'evals', out.funcCount);
if r.feasible, r.J = z(1); end
end

function y = logit(p), p = min(max(p, 1e-4), 1-1e-4); y = log(p/(1-p)); end
