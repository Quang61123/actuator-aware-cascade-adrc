function [A, S] = lockstep_miss_sensitivity()
%LOCKSTEP_MISS_SENSITIVITY Lockstep misses of the redesign SITL campaign.
%   (a) Each run's wall-clock interval [end - wall, end] (end = save time of
%       its result file, wall = run.wall) is intersected with the intervals
%       in which another MATLAB job ran on the machine during the campaign
%       (declared from the job logs, before looking at the outcome):
%         J1 18:07:30-18:11:40  design_pid_budget_curve, 11-worker pool (stopped by hand)
%         J2 18:15:40-18:18:30  est_design_iteration(0) refit, 1 worker
%       on 14-Oct-2026.
%   (b) Paired comparisons of the -E designs (S6, S5) with all runs and with
%       every pair that contains a run with a miss removed.
%   (c) Decision rule, fixed here: if a run of the main comparison groups
%       (EFS: ids 17-19 in S6) lies in J1 or J2, the whole group within that
%       interval is rerun and both sets are reported; otherwise no rerun.
%   Writes ../results/analysis/miss_by_run.csv, miss_by_group.csv,
%   lockstep_miss_sensitivity.csv.
here = fileparts(mfilename('fullpath')); RA = fullfile(here,'..','results','analysis'); RR = fullfile(RA,'runs');
J = [datetime(2026,10,14,18,7,30) datetime(2026,10,14,18,11,40); datetime(2026,10,14,18,15,40) datetime(2026,10,14,18,18,30)];
f = dir(fullfile(RR,'*.mat')); A = table();
for k = 1:numel(f)
    Z = load(fullfile(RR, f(k).name)); if ~isfield(Z,'run'), continue, end
    t1 = datetime(f(k).datenum, 'ConvertFrom', 'datenum'); t0 = t1 - seconds(Z.run.wall);
    ov = any(t0 < J(:,2) & t1 > J(:,1));
    A = [A; table(string(Z.run.group), Z.run.ctrl, Z.run.seed, Z.M.rmse3*1e3, Z.X.J_unq, Z.run.miss, t0, t1, ov, ...
        'VariableNames', {'Group','Ctrl','Seed','RMSE_mm','J_unq','Miss','Start','End','InParallelJob'})]; %#ok<AGROW>
end
writetable(A, fullfile(RA,'miss_by_run.csv'));
G = groupsummary(A, {'Group','Ctrl'}, {'sum','max'}, {'Miss'});
G.RunsWithMiss = arrayfun(@(i) sum(A.Group == G.Group(i) & A.Ctrl == G.Ctrl(i) & A.Miss > 0), (1:height(G))');
G.RunsInParallelJob = arrayfun(@(i) sum(A.Group == G.Group(i) & A.Ctrl == G.Ctrl(i) & A.InParallelJob), (1:height(G))');
writetable(G, fullfile(RA,'miss_by_group.csv'));
% (b) sensitivity of the paired -E comparisons
REST = fullfile(RA,'..','estimator_noise','runs'); S = table();
for sc = ["EFS","EF5"]
    if sc == "EFS", seeds = 91001:91030; pp = 'est_E_s6_c1_seed%d.mat'; else, seeds = 91001:91010; pp = 'est_E_s5_c1_seed%d.mat'; end
    P = nan(numel(seeds),2); for q = 1:numel(seeds), Z = load(fullfile(REST, sprintf(pp, seeds(q)))); P(q,:) = [Z.M.rmse3*1e3 Z.run.miss]; end
    get = @(id) arrayfun(@(s) pickv(A, sc, id, s), seeds, 'uni', 0);
    X = cell(1,19); for id = 17:19, c = get(id); X{id} = cell2mat(c(:)); end
    for id = [18 19]
        for ref = ["PID-E","PID-Budget"]
            if ref == "PID-E", y = X{17}(:,1); my = X{17}(:,2); else, y = P(:,1); my = P(:,2); end
            x = X{id}(:,1); mx = X{id}(:,2);
            for mode = ["all","no-miss pairs"]
                k = true(numel(seeds),1); if mode == "no-miss pairs", k = mx == 0 & my == 0; end
                d = x(k) - y(k); n = numel(d); t = tinv(0.975, max(n-1,1));
                S = [S; table(sc, id, ref, mode, n, mean(d), t*std(d)/sqrt(n), sum(d < 0), ...
                    'VariableNames', {'Scenario','Ctrl','Reference','Pairs','N','MeanDiff_mm','CI_mm','Better'})]; %#ok<AGROW>
            end
        end
    end
end
writetable(S, fullfile(RA,'lockstep_miss_sensitivity.csv'));
end

function v = pickv(A, g, id, s)
r = A(A.Group == g & A.Ctrl == id & A.Seed == s, :);
if isempty(r), v = [NaN NaN]; else, v = [r.RMSE_mm(1) r.Miss(1)]; end
end
