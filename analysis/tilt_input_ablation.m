function T = tilt_input_ablation()
%TILT_INPUT_ABLATION Ablation of the horizontal observer input: AA-Budget whose x-y ESOs receive the
%   commanded tilt (g*theta_d, g*phi_d after clipping) instead of the measured attitude (kind 10; inner and
%   altitude classes unchanged, actuator-aware). (i) fixed AA-Budget gains; (ii) horizontal class reselected
%   with the procedure of AA-Budget (cascade_outer_design: 20 ratios 1.5-60, local ceilings, Ms, OS,
%   primary cascade ceiling), then audited with verify_design (incl. settling). MATLAB S6 (91001-91030) and
%   S5 (91001-91010), non-MEX simulator (kind 10 is not in the compiled MEX). Paired against AA-Budget
%   (trials). Writes tilt_verify.csv, tilt_runs.csv, tilt_summary.csv, tilt_designs.mat.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); P = sim_config(); D = load(fullfile(RP,'design_gains.mat')); S = D.S;
z0 = cascade_stats(D.C(2), S); fprintf('regression check AA-Budget cascade du %.6f (reference 0.937)\n', z0(4));
A = D.C(2); A.kind = 10; A.name = 'A4-CommandedTilt (fixed gains)'; DS = A; V = verify_design(A, S, D.dub, [15 5 2 5], D.cap);
try
  [B, led] = cascade_outer_design(A, S, D.dub, D.cap); B.name = 'A4-CommandedTilt (retuned)'; DS = [DS, B]; V = [V; verify_design(B, S, D.dub, [15 5 2 5], D.cap)]; writetable(led, fullfile(RA,'tilt_ledger.csv'));
catch e, fprintf('RETUNE_FAILED %s\n', e.message); end
writetable(V, fullfile(RA,'tilt_verify.csv')); save(fullfile(RA,'tilt_designs.mat'), 'DS'); disp(V)
R = table(); cases = {6, 30, 40, [10 40]; 5, 10, 45, [5 45]};
for d = 1:numel(DS), co = DS(d); for s = 1:2, sc = cases{s,1}; n = cases{s,2}; Te = cases{s,3}; win = cases{s,4}; Z = nan(n,5);
  parfor q = 1:n
    seed = P.nominalSeed + q; pv = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat); E = make_exo(sc, seed, Te, S.Ts, S);
    L = run_quad_stress(co.kind, co.G, sc, pv, pc, Te, S.Ts, E, S.sig, 5*ones(round(Te/S.Ts),1)); m = run_metrics(L, win, pv);
    Z(q,:) = [seed, m.rmse3*1000, m.max3*1000, m.chat_cmd, m.diverged];
  end
  t = array2table(Z, 'VariableNames', {'Seed','RMSE_mm','MaxError_mm','Increment','Failed'}); t.Design = repmat(string(co.name), n, 1); t.Scen = repmat(sc, n, 1); R = [R; t]; %#ok<AGROW>
end, end
writetable(R, fullfile(RA,'tilt_runs.csv'));
O = readtable(fullfile(RP,'trials.csv'),'TextType','string'); O = O(abs(O.Tau-0.03)<1e-9 & O.Delay==5 & O.Noise==1 & O.Pattern=="fixed" & O.Ctrl=="AA-Budget",:); T = table();
for d = 1:numel(DS), for sc = [6 5]
  a = R(R.Design==string(DS(d).name) & R.Scen==sc,:); r = O(O.Scenario==sc,:); [~,ia] = unique(r.Seed,'first'); r = r(ia,:); [~,i1,i2] = intersect(a.Seed, r.Seed);
  x = a.RMSE_mm(i1) - r.RMSE_mm(i2); n = numel(x);
  T = [T; table(string(DS(d).name), sc, n, mean(a.RMSE_mm), std(a.RMSE_mm), mean(r.RMSE_mm(i2)), mean(x), tinv(0.975,n-1)*std(x)/sqrt(n), sum(x>0), mean(a.Increment), sum(a.Failed), 'VariableNames', {'Design','Scen','N','RMSE','SD','AAB_RMSE','MinusAAB','CI','AABBetterPairs','Inc','Failed'})]; %#ok<AGROW>
end, end
writetable(T, fullfile(RA,'tilt_summary.csv')); disp(T)
end