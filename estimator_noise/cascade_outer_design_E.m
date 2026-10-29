function [co, ledger] = cascade_outer_design_E(co, S, dub, capE, NF)
%CASCADE_OUTER_DESIGN_E Outer x-y redesign under an estimator-aware cascade
%   ceiling. Same search as cascade_outer_design (observer
%   ratios 1.5-60, bandwidth grid and bisection, cost = tracking standard
%   deviation, full step check of the selected candidates) with one change:
%   the cascade-level inner-increment ceiling is evaluated with the
%   estimator-coloured noise model (est_stats) against capE, and the cost
%   is the estimator-aware tracking standard deviation. The local outer
%   ceilings (white model, as in the primary design), Ms <= S.Msmax and OS <= S.osmax
%   are unchanged. Returns co = [] if no candidate is feasible.
original = co; ratios = logspace(log10(1.5), log10(60), 20); if co.kind == 2, ratios = 1; end
candidates = cell(numel(ratios)+1,1); cost = inf(numel(ratios)+1,1); ledger = table();
candidates{end} = original; if valid(original, true), z = est_stats(original, S, NF); cost(end) = z(1); end
for j = 1:numel(ratios)
    cc = original; shape = ratios(j); upper = S.wo_max/shape;
    grid = logspace(log10(.05), log10(upper), 36); best = NaN; idx = 0;
    for i = numel(grid):-1:1
        cc.G(1:6) = gains_of(cc.kind, shape, grid(i))'; cc.G(7:12) = cc.G(1:6);
        if valid(cc, false), best = grid(i); idx = i; break, end
    end
    if isfinite(best)
        if idx < numel(grid)
            lo = log(best); hi = log(grid(idx+1));
            for b = 1:14
                w = exp((lo+hi)/2); cc.G(1:6) = gains_of(cc.kind, shape, w)'; cc.G(7:12) = cc.G(1:6);
                if valid(cc, false), lo = log(w); else, hi = log(w); end
            end
            best = exp(lo);
        end
        cc.G(1:6) = gains_of(cc.kind, shape, best)'; cc.G(7:12) = cc.G(1:6);
        if valid(cc, true), zz = est_stats(cc, S, NF); cost(j) = zz(1); candidates{j} = cc; end
    end
end
[m, j] = min(cost);
if ~isfinite(m), co = []; return, end
co = candidates{j};
for j = 1:numel(candidates)
    if ~isempty(candidates{j})
        zz = est_stats(candidates{j}, S, NF);
        ledger = [ledger; table(string(co.name), capE, j, zz(1), zz(4), isfinite(cost(j)), ...
            'VariableNames', {'Ctrl','CascadeCapE','Candidate','SigmaE_est','CascadeSigmaDU_est','Feasible'})]; %#ok<AGROW>
    end
end
    function ok = valid(c, full)
        zw = cascade_stats(c, S);
        ok = all(isfinite(zw)) && zw(2) <= S.ub(1)*(1+1e-9) && zw(3) <= dub(1)*(1+1e-9);
        if ~ok, return, end
        ze = est_stats(c, S, NF); ok = all(isfinite(ze)) && ze(4) <= capE*(1-0.01);
        if ~ok, return, end
        ms = chan_ms(c.kind, c.G(1:6)', S.Ts, c.G(19:24)', S.act, [], [c.tauhat c.ndhat], 600); ok = ms <= S.Msmax;
        if ~ok, return, end
        if full, os = chan_step(c.kind, c.G(1:6)', S.Ts, c.G(19:24)', S.act, 4*S.Tstep(1), [c.tauhat c.ndhat]);
        else, os = step_screen(c.kind, c.G(1:6)', S.Ts, c.G(19:24)', S.act, 4*S.Tstep(1), [c.tauhat c.ndhat]); end
        ok = os <= S.osmax;
    end
end
