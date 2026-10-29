function T = make_run_plan_setpoint()
%MAKE_RUN_PLAN_SETPOINT Fixed SITL plan of phase P4: PID-Budget (b*), controller id 22.
%   PBS: S6 seeds 91001-91030 (EKF2 input); PB5: S5 seeds 91001-91010; PBR: S6 seeds 91001-91010 with the
%   estimator bypassed (raw measurement, as the act_raw diagnostic). RunIDs 2001-2050. Writes run_plan_setpoint.csv.
here = fileparts(mfilename('fullpath')); RA = fullfile(here,'..','results','analysis'); T = table();
G = {"PBS","S6",91001:91030,40,[10 40]; "PB5","S5",91001:91010,45,[5 45]; "PBR","S6",91001:91010,40,[10 40]}; id = 2000;
for g = 1:3, for s = G{g,3}, id = id + 1; w = G{g,5};
  T = [T; table(id, "P4", G{g,1}, "22", G{g,2}, s, G{g,4}, w(1), w(2), NaN, 'VariableNames', {'RunID','Phase','Group','CtrlId','Scenario','Seed','Tend_s','WinStart_s','WinEnd_s','Split'})]; %#ok<AGROW>
end, end
writetable(T, fullfile(RA,'run_plan_setpoint.csv')); disp(height(T))
end