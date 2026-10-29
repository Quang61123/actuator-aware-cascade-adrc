function T = sitl_collect()
%SITL_COLLECT Gather all saved SITL runs into results/sitl_all_runs.csv.
here = fileparts(mfilename('fullpath'));
R = fullfile(here,'results','runs'); d = dir(fullfile(R,'*.mat'));
G = load(fullfile(here,'sitl_gain_table.mat'),'names'); names = G.names;
T = table();
for i = 1:numel(d)
    S = load(fullfile(R,d(i).name),'M','run'); M = S.M; r = S.run;
    T = [T; table(string(names{r.ctrl}), r.ctrl, r.scen, r.seed, r.tau_ms, r.rep, r.diverged, r.miss, r.wall, ...
        M.rmse3*1e3, M.max3*1e3, M.rmsez*1e3, M.chat_cmd, M.rotor_std, M.sat_rotor, ...
        'VariableNames',{'Ctrl','id','Scen','Seed','TauPlant_ms','Rep','Failed','Miss','Wall_s', ...
        'RMSE_mm','Max_mm','RMSEz_mm','Increment','RotorStd','SatRotor'})]; %#ok<AGROW>
end
if ~isempty(T)
    T.Threshold = T.RotorStd > 20 & ~T.Failed;
    T = sortrows(T, {'Scen','TauPlant_ms','Seed','Rep','id'});
    writetable(T, fullfile(here,'results','sitl_all_runs.csv'));
end
end
