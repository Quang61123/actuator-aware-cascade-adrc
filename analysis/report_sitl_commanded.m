function T = report_sitl_commanded()
%REPORT_SITL_COMMANDED Phase P3 (run_plan_commanded.csv) SITL results and paired SITL comparisons (sitl_vs_matlab table).
%   A1X (id 20, A1-Commanded ratio>=1.0) S6 91001-010 vs AA-Budget; ROS (id 21, RO-HalfGain ratio>=1.0)
%   S6 91001-030 vs AA-HalfGain; RO5 (id 21) S5 91001-010 vs AA-HalfGain; PSS (id 15, new PID-Budget-S)
%   S6 91001-030 vs AA-Budget-S / AA-HalfGain-S. Reference runs: PX4_SITL/results/sitl_all_runs.csv.
%   All runs reported; writes sitl_P3.csv.
here = fileparts(mfilename('fullpath')); RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
A = readtable(fullfile(RA,'sitl_all_runs.csv'),'TextType','string');
O = readtable(fullfile(here,'..','PX4_SITL','results','sitl_all_runs.csv'),'TextType','string'); O = O(O.TauPlant_ms==0 & O.Rep==0,:);
C = {"A1X",20,"A1-Commanded (ratio>=1.0)",6,91001:91010,"AA-Budget"; "ROS",21,"RO-HalfGain (ratio>=1.0)",6,91001:91030,"AA-HalfGain"; "RO5",21,"RO-HalfGain (ratio>=1.0)",5,91001:91010,"AA-HalfGain"; "PSS",15,"PID-Budget-S (new)",6,91001:91030,"AA-Budget-Strict"; "PSS",15,"PID-Budget-S (new)",6,91001:91030,"AA-HalfGain-Strict"};
T = table();
for c = 1:size(C,1)
  s = C{c,5}; x = nan(numel(s),1); ms = zeros(numel(s),1); fl = false(numel(s),1); r = nan(numel(s),1);
  B = A(A.Group==C{c,1} & A.Ctrl==C{c,2},:); [~,ia,ib] = intersect(s,B.Seed); x(ia) = B.RMSE_mm(ib); ms(ia) = B.Miss(ib); fl(ia) = B.Failed(ib);
  R = O(O.Ctrl==C{c,6} & O.Scen==C{c,4},:); [~,ia,ib] = intersect(s,R.Seed); r(ia) = R.RMSE_mm(ib);
  d = r - x; ok = ~isnan(d) & ~fl; n = sum(ok); t = tinv(0.975,n-1); okm = ok & ms==0; nm = sum(okm);
  T = [T; table(string(C{c,1}), string(C{c,3}), C{c,4}, sum(~isnan(x)), mean(x,'omitnan'), std(x,'omitnan'), sum(fl), sum(ms), sum(ms>0), string(C{c,6}), mean(r(ok)), mean(d(ok)), t*std(d(ok))/sqrt(n), sum(d(ok)<0), nm, mean(d(okm)), tinv(0.975,max(nm-1,1))*std(d(okm))/sqrt(nm), ...
     'VariableNames',{'Group','Design','Scen','N','RMSE_mm','RMSE_sd','Failed','MissedSteps','RunsWithMiss','Reference','Ref_RMSE_mm','RefMinusDesign_mm','CI_mm','RefBetterPairs','NoMissPairs','RefMinusDesign_noMiss','CI_noMiss'})]; %#ok<AGROW>
end
writetable(T, fullfile(RA,'sitl_P3.csv')); disp(T)
end