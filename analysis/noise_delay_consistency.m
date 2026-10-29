function [V, T] = noise_delay_consistency()
%NOISE_DELAY_CONSISTENCY Design-model consistency for the yaw rate. The -E designs were selected with the
%   estimator-campaign pooled noise model, whose yaw-rate errors were defined with a
%   0-ms delay, while the stack delay convention uses 6 ms for all body rates.
%   Here the pooled model is refitted with the stack convention on exactly the
%   data of the estimator-campaign fit (hover logs of its six designs, seeds
%   96001-96005 for fitting), i.e. without any log of the -E designs, and
%   (i) all constraints of the -E designs are re-checked with it
%       (verify_design_stack, Ms on 4000 frequencies);
%   (ii) the out-of-sample predictions (level c) of the -E designs are
%       recomputed against the 300-s hover runs (fit and held-out seed sets).
%   Writes ../results/analysis/consistent_verify.csv,
%   consistent_levelc.csv and noise_est_data_rateconv/.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); F = load(fullfile(RA,'frozen_ceilings_and_delays.mat')); ds = F.ds;
E7 = load(fullfile(RA,'designs_E.mat'));
nd = fullfile(RA,'noise_est_data_rateconv');
if exist(fullfile(nd,'noise_models.mat'),'file'), Ln = load(fullfile(nd,'noise_models.mat')); NFc = Ln.NF;
else
    f4 = dir(fullfile(px4sitl_work_dir(),'est_data','est_hover_*.mat'));
    NFc = fit_noise_models(fullfile({f4.folder}, {f4.name}), ds, [20 60], 96001:96005, 96006:96010, nd);
end
V = table();
for j = 1:3, V = [V; verify_design_stack(E7.DE(j), D.S, D.dub, NFc, ds, F.Bcas, F.Bpsi, 4000)]; end %#ok<AGROW>
writetable(V, fullfile(RA,'consistent_verify.csv'));
% level (c) with the consistent model against the 300-s runs
Tv = readtable(fullfile(RA,'validation_3level.csv'), 'TextType', 'string');
T = table(); ids = [17 18 19];
names = ["se_x","se_y","se_z","se_psi","su_phi","su_theta","su_psi","su_z","sdu_phi","sdu_theta","sdu_psi","sdu_z","J"];
for j = 1:3
    z = predict_estimator_stats(E7.DE(j), NFc, ds, 'pooled');
    p = [z.se_x z.se_y z.se_z z.se_psi z.su_phi z.su_theta z.su_psi z.su_z z.sdu_phi z.sdu_theta z.sdu_psi z.sdu_z z.J];
    for q = 1:13
        r = Tv(Tv.CtrlId == ids(j) & Tv.Quantity == names(q), :);
        T = [T; table(string(E7.DE(j).name), names(q), r.Measured_fit, r.Measured_heldout, r.Pred_c_designtime, p(q), ...
            p(q)/r.Measured_fit, p(q)/r.Measured_heldout, 'VariableNames', ...
            {'Design','Quantity','Measured_fit','Measured_heldout','Pred_c_old','Pred_c_consistent','Rc_fit','Rc_heldout'})]; %#ok<AGROW>
    end
end
writetable(T, fullfile(RA,'consistent_levelc.csv'));
end
