function [T, best] = tilt_input_scan()
%TILT_INPUT_SCAN Retuning of A4-CommandedTilt under the FULL requirement set, including the
%   response requirement: the cascade_outer_design search used in tilt_input_ablation does not check settling, and
%   its best design settles in 15.4 s (> 15 s). Here the horizontal class of A4 (kind 10: commanded-tilt
%   x-y ESO input, AA-Budget inner/altitude gains fixed) is searched on the grid of the commanded-input
%   scans (ratios 1.0-1.36 and the 20 ratios 1.5-60, 60 controller bandwidths each) against sigma_u <= ub,
%   sigma_du <= dub (local), Ms <= 1.98, rho < 1, OS <= 10 %, ts <= 15 s and the primary cascade ceiling
%   (1 % margin); the feasible point with the smallest sigma_e is audited (verify_design) and run in MATLAB
%   S6 (91001-91030) and S5 (91001-91010) with the non-MEX simulator. Writes tilt_scan_grid.csv,
%   tilt_scan_verify.csv, tilt_scan_runs.csv, tilt_scan_summary.csv.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); P = sim_config();
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub; cap = D.cap;
K1 = load(fullfile(RA,'commanded_scan_ratio1.mat')); rg = K1.rg(:)'; nw = 60;
base = D.C(2); base.kind = 10; base.name = 'A4-CommandedTilt (retuned, full requirements)'; gin = base.G(19:24)'; oa = [base.tauhat base.ndhat];
rows = cell(numel(rg),1);
parfor q = 1:numel(rg)
    r = rg(q); wg = logspace(log10(0.05), log10(S.wo_max/r), nw); R = nan(nw, 11);
    for i = 1:nw
        g = gains_of(10, r, wg(i)); z = channel_stats(10, g, 1, S, gin, oa, S.act, false, 600);
        os = NaN; ts = NaN; cdu = NaN; ok = false;
        if all(isfinite(z(1:4))) && z(2) <= S.ub(1)*(1+1e-9) && z(3) <= dub(1)*(1+1e-9) && z(4) <= S.Msmax && z(6) < 1
            [os, ts] = step_response_metrics(chan_model(10, g, S.Ts, gin, S.act, [], oa), S.Ts);
            c2 = base; c2.G(1:6) = g'; c2.G(7:12) = g'; zz = cascade_stats(c2, S); cdu = zz(4);
            ok = os <= S.osmax && ts <= 15 && cdu <= cap*(1-0.01);
        end
        R(i,:) = [r, wg(i), r*wg(i), z(1), z(2), z(3), z(4), os, ts, cdu, ok];
    end
    rows{q} = R;
end
G = array2table(vertcat(rows{:}), 'VariableNames', {'ratio','wc','wo','sigma_e','sigma_u','sigma_du','Ms','OS','ts','cascade_du','Feasible'});
writetable(G, fullfile(RA,'tilt_scan_grid.csv')); fe = G(G.Feasible == 1, :);
fprintf('C3SCAN feasible %d of %d\n', height(fe), height(G));
T = table(); best = [];
if isempty(fe), fprintf('C3SCAN no A4 design satisfies all requirements on this grid\n'); return, end
[~, ord] = sort(fe.sigma_e);
for j = ord'   % fine-grid Ms re-check, as in commanded_input_scan
    g = gains_of(10, fe.ratio(j), fe.wc(j)); if chan_ms(10, g, S.Ts, gin, S.act, [], oa, 4000) <= S.Msmax, best = base; best.G(1:6) = g'; best.G(7:12) = g'; break, end
end
if isempty(best), fprintf('C3SCAN no candidate passes the 4000-point Ms check\n'); return, end
V = verify_design(best, S, dub, [15 5 2 5], cap); writetable(V, fullfile(RA,'tilt_scan_verify.csv')); disp(V)
co = best; R = table(); cases = {6, 30, 40, [10 40]; 5, 10, 45, [5 45]};
for s = 1:2, sc = cases{s,1}; n = cases{s,2}; Te = cases{s,3}; win = cases{s,4}; Z = nan(n,5);
    parfor q = 1:n
        seed = P.nominalSeed + q; pv = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat); E = make_exo(sc, seed, Te, S.Ts, S);
        L = run_quad_stress(co.kind, co.G, sc, pv, pc, Te, S.Ts, E, S.sig, 5*ones(round(Te/S.Ts),1)); m = run_metrics(L, win, pv);
        Z(q,:) = [seed, m.rmse3*1000, m.max3*1000, m.chat_cmd, m.diverged];
    end
    t = array2table(Z, 'VariableNames', {'Seed','RMSE_mm','MaxError_mm','Increment','Failed'}); t.Scen = repmat(sc, n, 1); R = [R; t]; %#ok<AGROW>
end
writetable(R, fullfile(RA,'tilt_scan_runs.csv'));
O = readtable(fullfile(RP,'trials.csv'),'TextType','string'); O = O(abs(O.Tau-0.03)<1e-9 & O.Delay==5 & O.Noise==1 & O.Pattern=="fixed" & O.Ctrl=="AA-Budget",:);
for sc = [6 5]
    a = R(R.Scen==sc,:); r = O(O.Scenario==sc,:); [~,ia] = unique(r.Seed,'first'); r = r(ia,:); [~,i1,i2] = intersect(a.Seed, r.Seed); x = a.RMSE_mm(i1) - r.RMSE_mm(i2); n = numel(x);
    T = [T; table(string(co.name), sc, n, mean(a.RMSE_mm), std(a.RMSE_mm), mean(r.RMSE_mm(i2)), mean(x), tinv(0.975,n-1)*std(x)/sqrt(n), sum(x>0), mean(a.Increment), sum(a.Failed), ...
        'VariableNames', {'Design','Scen','N','RMSE','SD','AAB_RMSE','MinusAAB','CI','AABBetterPairs','Inc','Failed'})]; %#ok<AGROW>
end
writetable(T, fullfile(RA,'tilt_scan_summary.csv')); disp(T)
end
