function [TE, TS, TH] = report_sitl()
%REPORT_SITL Redesign SITL results (run_plan.csv), all runs reported.
%   TE: -E designs (ids 17-19, final) versus PID-Budget behind EKF2:
%       S6 seeds 91001-91030 (PID-Budget: estimator-campaign runs est_E_s6_c1, same seeds,
%       same setup), S5 seeds 91001-91010 (PID-Budget: est_E_s5_c1): RMSE
%       mean/SD, paired change with 95 % CI and pairs better, unquantized J,
%       paired J change and pairs with lower J, failures, missed lockstep steps.
%   TH: noise-only hover (E0H, 60 s): increments per class and J, measured.
%   TS: PID-Budget-S (new, id 15) S6 30 seeds versus the stored SITL runs of
%       AA-Budget-S / AA-HalfGain-S; A1-Commanded (new, id 16) S6
%       10 seeds versus AA-Budget.
%   Writes ../results/analysis/sitl_E.csv, sitl_hover_E.csv,
%   sitl_S_A1.csv, sitl_all_runs.csv.
here = fileparts(mfilename('fullpath')); RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); RR = fullfile(RA,'runs');
REST = fullfile(RP,'estimator_noise','runs');
f = dir(fullfile(RR,'*.mat')); A = table();
for k = 1:numel(f)
    Z = load(fullfile(RR, f(k).name));
    if ~isfield(Z,'run'), A = [A; table(string(f(k).name), NaN, "", NaN, NaN, NaN, NaN, true, NaN, 'VariableNames', {'Tag','Ctrl','Group','Seed','RMSE_mm','J_unq','J_q','Failed','Miss'})]; continue, end %#ok<AGROW>
    J = NaN; if isfield(Z.X,'J_unq'), J = Z.X.J_unq; end
    A = [A; table(string(Z.run.tag), Z.run.ctrl, string(Z.run.group), Z.run.seed, Z.M.rmse3*1e3, J, Z.X.J_q, logical(Z.M.diverged), Z.run.miss, ...
        'VariableNames', {'Tag','Ctrl','Group','Seed','RMSE_mm','J_unq','J_q','Failed','Miss'})]; %#ok<AGROW>
end
writetable(A, fullfile(RA,'sitl_all_runs.csv'));
TE = table();
for sc = ["S6","S5"]
    if sc == "S6", grp = "EFS"; pp = 'est_E_s6_c1_seed%d.mat'; seeds = 91001:91030; else, grp = "EF5"; pp = 'est_E_s5_c1_seed%d.mat'; seeds = 91001:91010; end
    P = nan(numel(seeds), 2);
    for q = 1:numel(seeds), Z = load(fullfile(REST, sprintf(pp, seeds(q)))); P(q,:) = [Z.M.rmse3*1e3, Z.X.J_unq]; end
    TE = [TE; row("PID-Budget (estimator campaign)", sc, P(:,1), P(:,2), P, zeros(numel(seeds),1), zeros(numel(seeds),1))]; %#ok<AGROW>
    for id = 17:19
        B = A(A.Group == grp & A.Ctrl == id, :); [~, ia, ib] = intersect(seeds, B.Seed); X = nan(numel(seeds), 2); X(ia,:) = [B.RMSE_mm(ib) B.J_unq(ib)];
        fl = false(numel(seeds),1); fl(ia) = B.Failed(ib); ms = zeros(numel(seeds),1); ms(ia) = B.Miss(ib);
        TE = [TE; row(nm(id), sc, X(:,1), X(:,2), P, fl, ms)]; %#ok<AGROW>
    end
end
writetable(TE, fullfile(RA,'sitl_E.csv'));
TH = table();
for id = [17 18 19]
    rr = dir(fullfile(RR, sprintf('E0H_hover_c%d_*.mat', id))); v = [];
    for q = 1:numel(rr), Z = load(fullfile(RR, rr(q).name)); v(end+1,:) = [Z.X.sd_du Z.X.J_unq].^2; end %#ok<AGROW>
    s = sqrt(mean(v,1)); TH = [TH; table(nm(id), size(v,1), max(s(1:2)), s(3), s(4), s(5), 'VariableNames', {'Design','N','SdDU_rp','SdDU_psi','SdDU_z','J'})]; %#ok<AGROW>
end
writetable(TH, fullfile(RA,'sitl_hover_E.csv'));
TS = table();
B = A(A.Group == "PSS", :); TS = [TS; table("PID-Budget-S (new)", "S6", height(B), mean(B.RMSE_mm), std(B.RMSE_mm), mean(B.J_unq), sum(B.Failed), sum(B.Miss), 'VariableNames', {'Design','Scenario','N','RMSE_mm','RMSE_sd','J_unq','Failed','MissedSteps'})];
B = A(A.Group == "A1S", :); TS = [TS; table("A1-Commanded (new)", "S6", height(B), mean(B.RMSE_mm), std(B.RMSE_mm), mean(B.J_unq), sum(B.Failed), sum(B.Miss), 'VariableNames', TS.Properties.VariableNames)];
B = A(A.Group == "PSH", :); TS = [TS; table("PID-Budget-S (new) hover", "hover", height(B), mean(B.RMSE_mm), std(B.RMSE_mm), mean(B.J_unq), sum(B.Failed), sum(B.Miss), 'VariableNames', TS.Properties.VariableNames)];
writetable(TS, fullfile(RA,'sitl_S_A1.csv'));
end

function s = nm(id)
n = ["PID-E","AA-Budget-E","AA-HalfGain-E"]; s = n(id-16);
end

function t = row(name, sc, x, j, P, fl, ms)
d = x - P(:,1); dj = j - P(:,2); ok = ~isnan(d); n = sum(ok); tt = tinv(0.975, max(n-1,1));
t = table(string(name), sc, n, mean(x,'omitnan'), std(x,'omitnan'), mean(d(ok)), tt*std(d(ok))/sqrt(n), sum(d(ok) < 0), ...
    mean(j,'omitnan'), mean(dj(ok)), tt*std(dj(ok))/sqrt(n), sum(dj(ok) < 0), sum(fl), sum(ms), ...
    'VariableNames', {'Design','Scenario','N','RMSE_mm','RMSE_sd','DiffVsPID_mm','DiffCI_mm','PairsBetter','J_unq','DiffJ','DiffJ_CI','PairsLowerJ','Failed','MissedSteps'});
end
