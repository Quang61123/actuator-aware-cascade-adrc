function T = make_run_plan_commanded()
%MAKE_RUN_PLAN_COMMANDED Fixed list of the additional redesign SITL runs (written
%   before running; run_plan_commanded.csv), needed because the commanded-input
%   designs changed after the ratio range was extended to 1.0:
%     A1X  S6 40 s, id 20 A1-Commanded (ratio >= 1.0), seeds 91001-91010
%     ROS  S6 40 s, id 21 RO-HalfGain (ratio >= 1.0), seeds 91001-91030
%     RO5  S5 45 s, id 21, seeds 91001-91010
%   Same runner (run_sitl_redesign, phase "P3"); every listed run is executed once
%   and reported; no other job runs on the machine during the campaign.
here = fileparts(mfilename('fullpath')); RA = fullfile(here,'..','results','analysis');
rows = cell(0, 9);
for s = 91001:91010, rows(end+1,:) = {"P3","A1X","20","S6",s,40,10,40,""}; end %#ok<AGROW>
for s = 91001:91030, rows(end+1,:) = {"P3","ROS","21","S6",s,40,10,40,""}; end %#ok<AGROW>
for s = 91001:91010, rows(end+1,:) = {"P3","RO5","21","S5",s,45,5,45,""}; end %#ok<AGROW>
T = cell2table(rows, 'VariableNames', {'Phase','Group','CtrlId','Scenario','Seed','Tend_s','WinStart_s','WinEnd_s','Split'});
T.RunID = 1000 + (1:height(T))'; T = T(:, [end 1:end-1]);
writetable(T, fullfile(RA,'run_plan_commanded.csv'));
end
