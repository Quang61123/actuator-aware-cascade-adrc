function T = control_std_check()
%CONTROL_STD_CHECK Measured control standard deviations of the -E designs against
%   the local variance ceilings ub(c) of the covariance design problem.
%   The local ceilings of the -E search are evaluated with the white-noise
%   design model (conservative for the increments, which it overpredicts
%   19-410 fold), whereas the estimator-aware model underpredicts sigma_u of
%   roll/pitch by up to 19 % (non-conservative). Here the logged sigma_u of
%   [u_phi u_theta u_psi u_z] in noise-only hover (300 s, seeds 97001-97010,
%   window 20-300 s; 60 s, seeds 96001-96010) is compared with ub(c), and the
%   white-model and estimator-aware predictions are listed. The horizontal
%   virtual control is not logged and is reported by prediction only.
%   Writes ../results/analysis/control_std_check.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); RR = fullfile(RA,'runs');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; F = load(fullfile(RA,'frozen_ceilings_and_delays.mat'));
E7 = load(fullfile(RA,'designs_E.mat')); Ln = load(fullfile(RA,'noise_long','noise_models.mat'));
ids = [1 17 18 19]; des = {D.C(1), E7.DE(1), E7.DE(2), E7.DE(3)}; T = table();
cls = ["phi","theta","psi","z"]; ci = [3 3 4 2];       % class index of ub for [u_phi u_theta u_psi u_z]
for j = 1:4
    co = des{j}; id = ids(j);
    rl = [dir(fullfile(RR, sprintf('*_hover_c%d_seed97*_T300.mat', id)))]; rs = dir(fullfile(RR, sprintf('E0H_hover_c%d_seed*_T60.mat', id)));
    if id == 1, rs = dir(fullfile(RP,'estimator_noise','runs','est_hover_c1_seed*.mat')); end
    vl = []; for q = 1:numel(rl), Z = load(fullfile(rl(q).folder, rl(q).name)); vl(end+1,:) = Z.X.sd_u.^2; end %#ok<AGROW>
    vs = []; for q = 1:numel(rs), Z = load(fullfile(rs(q).folder, rs(q).name)); vs(end+1,:) = Z.X.sd_u.^2; end %#ok<AGROW>
    ml = sqrt(mean(vl,1)); ms = sqrt(mean(vs,1));
    zb = predict_estimator_stats(co, Ln.NF, F.ds, 'pooled'); pe = [zb.su_phi zb.su_theta zb.su_psi zb.su_z];
    gs = {co.G(19:24)', co.G(19:24)', co.G(31:36)', co.G(13:18)'}; pw = zeros(1,4);
    for k = 1:4, q = channel_stats(co.kind, gs{k}, ci(k), S, [], [co.tauhat co.ndhat], S.act, false, 600); pw(k) = q(2); end
    for k = 1:4
        T = [T; table(string(co.name), cls(k), S.ub(ci(k)), ml(k), ms(k), pw(k), pe(k), max(ml(k), ms(k)) <= S.ub(ci(k)), ...
            'VariableNames', {'Design','Control','Ceiling_ub','Meas_hover300','Meas_hover60','Pred_white_local','Pred_est_pooled_cascade','MeasuredWithinCeiling'})]; %#ok<AGROW>
    end
    q = channel_stats(co.kind, co.G(1:6)', 1, S, co.G(19:24)', [co.tauhat co.ndhat], S.act, false, 600); zc = cascade_stats(co, S);
    T = [T; table(string(co.name), "xy (outer)", S.ub(1), NaN, NaN, zc(2), NaN, zc(2) <= S.ub(1)*(1+1e-9), 'VariableNames', T.Properties.VariableNames)]; %#ok<AGROW>
end
writetable(T, fullfile(RA,'control_std_check.csv'));
end
