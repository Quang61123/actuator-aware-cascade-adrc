function T = commanded_low_ratio()
%COMMANDED_LOW_RATIO Sensitivity of the commanded-input designs to the observer/controller ratio
%   range: A1-Commanded and RO-HalfGain searched on the ratios
%   omega_o/omega_c in {0.6, 0.8} only (all constraints of the covariance design problem plus the response requirement), then MATLAB S6
%   (91001-91030) and paired against AA-Budget / AA-HalfGain (trials). Writes commanded_scan_ratio068.*,
%   lowratio_verify.csv, lowratio_runs.csv, lowratio_summary.csv.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); P = sim_config(); D = load(fullfile(RP,'design_gains.mat')); S = D.S;
rg = [0.6 0.8];   % designs restricted to omega_o/omega_c in {0.6, 0.8}
commanded_input_scan([1 2], rg, '_ratio068'); K = load(fullfile(RA,'commanded_scan_ratio068.mat'));
fam = {"A1-Commanded", D.C(5), "AA-Budget"; "RO-HalfGain", D.C(3), "AA-HalfGain"}; V = table(); R = table(); DS = {};
for f = 1:2
  A = K.T(string(K.T.Design) == fam{f,1}, :); co = fam{f,2}; cells = cell(1,4); ok = true;
  for c = 1:4, r = A(A.Class == c, :); if isempty(r) || ~isfinite(r.Best_wc), ok = false; break, end, cells{c} = [r.Best_wc r.Best_wo 0 0 0 0]; end
  if ~ok, fprintf('low-ratio %s: no feasible design on the extended grid\n', fam{f,1}); continue, end
  co.G = pack_gains(cells); co.name = char(fam{f,1} + " (ratio 0.6-0.8)"); DS{end+1} = co; V = [V; verify_design(co, S, D.dub, [15 5 2 5], D.cap)]; %#ok<AGROW>
  Z = nan(30,4);
  parfor q = 1:30
    seed = P.nominalSeed + q; pv = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat); E = make_exo(6, seed, 40, S.Ts, S);
    L = run_quad_stress_mex(co.kind, co.G, 6, pv, pc, 40, S.Ts, E, S.sig, 5*ones(40000,1)); m = run_metrics(L, [10 40], pv); Z(q,:) = [seed, m.rmse3*1000, m.chat_cmd, m.diverged];
  end
  t = array2table(Z, 'VariableNames', {'Seed','RMSE_mm','Increment','Failed'}); t.Design = repmat(string(co.name),30,1); t.Ref = repmat(fam{f,3},30,1); R = [R; t]; %#ok<AGROW>
end
writetable(V, fullfile(RA,'lowratio_verify.csv')); writetable(R, fullfile(RA,'lowratio_runs.csv')); disp(V)
O = readtable(fullfile(RP,'trials.csv'),'TextType','string'); O = O(O.Experiment=="E1",:); T = table();
for d = unique(R.Design)'
  a = R(R.Design==d,:); ref = a.Ref(1); r = O(O.Ctrl==ref,:); [~,i1,i2] = intersect(a.Seed, r.Seed); x = a.RMSE_mm(i1) - r.RMSE_mm(i2); n = numel(x);
  T = [T; table(d, ref, n, mean(a.RMSE_mm), std(a.RMSE_mm), mean(r.RMSE_mm(i2)), mean(x), tinv(0.975,n-1)*std(x)/sqrt(n), sum(x>0), mean(a.RMSE_mm)/mean(r.RMSE_mm(i2)), mean(a.Increment), 'VariableNames', {'Design','Ref','N','RMSE','SD','RefRMSE','Diff','CI','RefBetter','Ratio','Inc'})]; %#ok<AGROW>
end
writetable(T, fullfile(RA,'lowratio_summary.csv')); disp(T)
end