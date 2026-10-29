function T = matlab_paired()
%MATLAB_PAIRED Paired MATLAB comparisons of the redesigned designs with the reference designs.
%   Reference per-seed results: results/trials.csv (nominal plant: tau 30 ms, delay 5,
%   noise 1, fixed pattern; seeds 91001-91030). New designs: matlab_<which>.csv (S6) and
%   matlab_ROx_S5.csv (S5). Difference = new - reference; CI = t-based 95% half-width. Writes matlab_paired.csv.
here = fileparts(mfilename('fullpath')); RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
O = readtable(fullfile(RP,'trials.csv'),'TextType','string');
O = O(abs(O.Tau-0.03)<1e-9 & O.Delay==5 & O.Noise==1 & O.Pattern=="fixed",:);
C = {'A1x',6,"AA-Budget"; 'ROx',6,"AA-HalfGain"; 'ROx_S5',5,"AA-HalfGain"; 'PIDS',6,"AA-Budget-Strict"; 'PIDS',6,"AA-HalfGain-Strict"; 'CURVE050',6,"P2-g050"; 'CURVE050',6,"HG-AA-g050"; 'CURVE075',6,"P2-g075"; 'CURVE075',6,"HG-AA-g075"};
T = table();
for c = 1:size(C,1)
  N = readtable(fullfile(RA,['matlab_' C{c,1} '.csv'])); R = O(O.Ctrl==C{c,3} & O.Scenario==C{c,2},:);
  [~,ia] = unique(R.Seed,'first'); R = R(ia,:); [s,ia,ib] = intersect(N.Seed, R.Seed);
  ok = ~N.Failed(ia) & ~R.Failed(ib); d = N.RMSE_mm(ia(ok)) - R.RMSE_mm(ib(ok)); n = numel(d);
  T = [T; table(string(C{c,1}), C{c,2}, string(C{c,3}), n, mean(N.RMSE_mm(ia(ok))), mean(R.RMSE_mm(ib(ok))), mean(R.Increment(ib(ok))), mean(d), tinv(0.975,n-1)*std(d)/sqrt(n), sum(d>0), mean(N.RMSE_mm(ia(ok)))/mean(R.RMSE_mm(ib(ok))), ...
    'VariableNames',{'New','Scen','Reference','N','New_RMSE','Ref_RMSE','Ref_Inc','NewMinusRef_mm','CI_mm','RefBetterPairs','RMSE_ratio'})]; %#ok<AGROW>
end
writetable(T, fullfile(RA,'matlab_paired.csv')); disp(T)
end