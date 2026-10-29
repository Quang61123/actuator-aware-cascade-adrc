function [T, pass, NFu] = est_design_iteration(iter)
%EST_DESIGN_ITERATION Redesign-log-refit check of the estimator-stack designs,
%   iteration iter (0: ids 17-19 from designs_E.mat,
%   1: ids 20-22 from designs_E_it1.mat; not used, since iteration 0 met the stopping rule).
%   1. Refit (fit_noise_models): hover logs (window 20-60 s, seeds 96001-96010,
%      fit on 96001-96005) of the six designs of the estimator campaign (est_data) and of the
%      -E designs of this iteration (redesign_data, group E<iter>H), with the stack
%      delay convention; per-design models and a model pooled over all
%      designs (the updated design model).
%   2. Predictions of sigma_du (roll/pitch, yaw, z) and J of the -E designs
%      with (c) the pooled model used at design time (estimator campaign, before any
%      log of these designs: out-of-sample), (b) the updated pooled model,
%      (a) the per-design models; measured values from the ten hover runs.
%   3. Stopping rule: |updated-pooled prediction / measured - 1| <= 10 % for
%      the roll/pitch and yaw increments of every -E design, and all
%      constraints of verify_design_stack still met with the updated pooled model.
%   Writes ../results/analysis/est_design_iteration<iter>.csv and noise_it<iter>/.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S;
F = load(fullfile(RA,'frozen_ceilings_and_delays.mat')); ds = F.ds;
L4 = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF0 = L4.NF;
if iter == 0, K = load(fullfile(RA,'designs_E.mat')); ids = 17:19; grp = 'E0H';
else, K = load(fullfile(RA,'designs_E_it1.mat')); ids = 20:22; grp = 'E1H'; end
DE = K.DE;
f4 = dir(fullfile(px4sitl_work_dir(),'est_data','est_hover_*.mat'));
f7 = dir(fullfile(fullfile(px4sitl_work_dir(),'redesign_data'), ['' grp '_hover_*.mat']));
files = [fullfile({f4.folder}, {f4.name}), fullfile({f7.folder}, {f7.name})];
nd = fullfile(RA, sprintf('noise_it%d', iter));
if exist(fullfile(nd,'noise_models.mat'),'file'), Ln = load(fullfile(nd,'noise_models.mat')); NFu = Ln.NF;   % reuse the refit of this iteration
else, NFu = fit_noise_models(files, ds, [20 60], 96001:96005, 96006:96010, nd); end
T = table(); pass = true;
for j = 1:3
    co = DE(j); id = ids(j);
    rr = dir(fullfile(RA,'runs', sprintf('%s_hover_c%d_seed*_T60.mat', grp, id))); m = [];
    for q = 1:numel(rr), Rq = load(fullfile(rr(q).folder, rr(q).name)); if isfield(Rq,'X') && isfield(Rq.X,'sd_du'), m(end+1,:) = [Rq.X.sd_du.^2, Rq.X.J_unq^2]; end, end %#ok<AGROW>
    meas = sqrt(mean(m,1));   % [phi theta psi z J]
    zc = predict_estimator_stats(co, NF0, ds, 'pooled'); zb = predict_estimator_stats(co, NFu, ds, 'pooled'); za = predict_estimator_stats(co, NFu, ds, 'design', id);
    q = @(z) [max(z.sdu_phi, z.sdu_theta), z.sdu_psi, z.sdu_z, z.J];
    mm = [max(meas(1:2)), meas(3), meas(4), meas(5)];
    V = verify_design_stack(co, S, D.dub, NFu, ds, F.Bcas, F.Bpsi, 4000);
    rb = q(zb)./mm;
    ok = all(abs(rb(1:2) - 1) <= 0.10) && all(V.Pass);
    pass = pass && ok;
    nm = {'sdu_rp','sdu_psi','sdu_z','J'};
    for k = 1:4
        T = [T; table(string(co.name), id, size(m,1), string(nm{k}), mm(k), q(zc)*[k==1;k==2;k==3;k==4], q(zb)*[k==1;k==2;k==3;k==4], q(za)*[k==1;k==2;k==3;k==4], ...
            'VariableNames', {'Design','CtrlId','Runs','Quantity','Measured','Pred_c_designtime','Pred_b_updatedPooled','Pred_a_perDesign'})]; %#ok<AGROW>
    end
    writetable(V, fullfile(RA, sprintf('est_design_iteration%d_verify_%d.csv', iter, id)));
    fprintf('ITER %d %s: rp %.3f psi %.3f constraints %d -> %d\n', iter, co.name, rb(1), rb(2), all(V.Pass), ok);
end
T.RatioB = T.Pred_b_updatedPooled ./ T.Measured; T.RatioC = T.Pred_c_designtime ./ T.Measured; T.RatioA = T.Pred_a_perDesign ./ T.Measured;
writetable(T, fullfile(RA, sprintf('est_design_iteration%d.csv', iter)));
save(fullfile(RA, sprintf('est_design_iteration%d.mat', iter)), 'T', 'pass');
end
