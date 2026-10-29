function [T6, TH, T5] = report_estimator_initial()
%REPORT_ESTIMATOR_INITIAL Equal-activity comparison in the estimator-fed loop.
%   PX4 SITL with EKF2 input, PID-Budget (unchanged) versus AA-Budget-E and
%   AA-HalfGain-E (design_estimator_initial, ids 11 and 12):
%     S6, seeds 91001-91030, window 10-40 s (est_E_s6_*): RMSE, unquantized
%       rotor-command increment J (debug_array), failures, paired RMSE
%       differences versus PID-Budget (mean, 95% CI half-width, pairs better)
%     noise-only hover, seeds 96001-96010 (est_E_hover_* and est_hover_c1_*):
%       inner increment sigma_du and J, measured versus the estimator-aware
%       prediction (pooled noise models)
%     S5, seeds 91001-91010, window 5-45 s (est_E_s5_*): RMSE, J
%   Writes est_E_s6.csv, est_E_hover.csv, est_E_s5.csv and est_E_s6_per_run.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RES = fullfile(RP,'estimator_noise'); RR = fullfile(RES,'runs');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; pv = quad_params();
L = load(fullfile(RES,'est_noise_models.mat')); NF = L.NF; E = load(fullfile(RES,'est_designs_E.mat'));
des = {1, D.C(1); 11, E.DE(1); 12, E.DE(2)};
% ---------------- S6 ----------------
A = table();
for r = 1:size(des,1)
    for q = 1:30
        f = fullfile(RR, sprintf('est_E_s6_c%d_seed%d.mat', des{r,1}, 91000+q)); if ~exist(f,'file'), continue, end
        R = load(f); A = [A; table(string(des{r,2}.name), des{r,1}, 91000+q, R.M.rmse3*1e3, R.M.max3*1e3, R.X.J_unq, R.X.J_q, ...
            double(R.M.diverged), R.run.miss, 'VariableNames', {'Design','Id','Seed','RMSE_mm','MaxErr_mm','J_unq','J_q','Diverged','Miss'})]; %#ok<AGROW>
    end
end
writetable(A, fullfile(RES,'est_E_s6_per_run.csv'));
T6 = table(); p = A(A.Id == 1,:);
for r = 1:size(des,1)
    a = A(A.Id == des{r,1},:); [~, ia, ip] = intersect(a.Seed, p.Seed); d = a.RMSE_mm(ia) - p.RMSE_mm(ip);
    ci = tinv(0.975, max(numel(d)-1,1))*std(d)/sqrt(numel(d));
    T6 = [T6; table(string(des{r,2}.name), height(a), mean(a.RMSE_mm), std(a.RMSE_mm), mean(a.J_unq), mean(a.J_q), sum(a.Diverged), ...
        mean(d), ci, sum(d < 0), 'VariableNames', {'Design','N','RMSE_mm','RMSE_sd','J_unq','J_q','Failed','DiffVsPID_mm','DiffCI_mm','PairsBetterThanPID'})]; %#ok<AGROW>
end
writetable(T6, fullfile(RES,'est_E_s6.csv'));
% ---------------- hover ----------------
TH = table();
for r = 1:size(des,1)
    id = des{r,1}; co = des{r,2};
    if id == 1, pat = 'est_hover_c1_seed%d.mat'; else, pat = ['est_E_hover_c' num2str(id) '_seed%d.mat']; end
    v = []; jj = [];
    for q = 1:10
        f = fullfile(RR, sprintf(pat, 96000+q)); if ~exist(f,'file'), continue, end
        R = load(f); v(end+1,:) = R.X.sd_du(1:2).^2; jj(end+1) = R.X.J_unq^2; %#ok<AGROW>
    end
    z = est_stats(co, S, NF); Jp = rotorJ_est(co, S, NF, pv);
    sm = max(sqrt(mean(v,1))); Jm = sqrt(mean(jj));
    TH = [TH; table(string(co.name), size(v,1), z(4), sm, sm/z(4), Jp, Jm, Jm/Jp, ...
        'VariableNames', {'Design','N','Pred_SigmaDU','SITL_SigmaDU','Ratio_SigmaDU','Pred_J','SITL_J','Ratio_J'})]; %#ok<AGROW>
end
writetable(TH, fullfile(RES,'est_E_hover.csv'));
% ---------------- S5 ----------------
T5 = table();
for r = 1:size(des,1)
    x = []; jj = []; dv = 0;
    for q = 1:10
        f = fullfile(RR, sprintf('est_E_s5_c%d_seed%d.mat', des{r,1}, 91000+q)); if ~exist(f,'file'), continue, end
        R = load(f); x(end+1) = R.M.rmse3*1e3; jj(end+1) = R.X.J_unq; dv = dv + double(R.M.diverged); %#ok<AGROW>
    end
    T5 = [T5; table(string(des{r,2}.name), numel(x), mean(x), std(x), mean(jj), dv, 'VariableNames', {'Design','N','RMSE_mm','RMSE_sd','J_unq','Failed'})]; %#ok<AGROW>
end
writetable(T5, fullfile(RES,'est_E_s5.csv'));
end

function J = rotorJ_est(co, S, NF, pv)
% estimator-aware rotor-command increment prediction (as predict_activity, pooled models)
ds = NF.dsig; oa = [co.tauhat co.ndhat]; Ts = S.Ts; act = S.act;
zX = cov_estimator_aware(chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, ds([1 4 8 11])), NF.X_pooled);
zY = cov_estimator_aware(chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, ds([2 5 7 10])), NF.Y_pooled);
zZ = cov_estimator_aware(chan_model_est(co.kind, co.G(13:18)', Ts, [], act, [], oa, [ds([3 6]) 0 0]), NF.Z_pooled);
zP = cov_estimator_aware(chan_model_est(co.kind, co.G(31:36)', Ts, [], act, [], oa, [ds([9 12]) 0 0]), NF.P_pooled);
k = pv(4); l = pv(3); b = pv(5); wh = sqrt(pv(1)*pv(2)/(4*k));
Mx = [k k k k; 0 -l*k 0 l*k; -l*k 0 l*k 0; -b b -b b]; Jm = diag(1/(2*wh)*ones(4,1))/Mx;
I3 = pv(7:9); v = [(pv(1)*zZ(3))^2, (I3(1)*zY(4))^2, (I3(2)*zX(4))^2, (I3(3)*zP(3))^2];
J = sqrt(mean((Jm.^2)*v(:)));
end
