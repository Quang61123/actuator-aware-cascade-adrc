function T = validate_noise_models()
%VALIDATE_NOISE_MODELS Three-level validation of the estimator-aware model for
%   every quantity used in the design, on 300-s
%   noise-only hover SITL runs (window 20-300 s; seeds 97001-97005 = fit,
%   97006-97010 = held-out, disjoint):
%     (a) per-design model fitted on the fit seeds of that design,
%     (b) model pooled over the nine designs, fitted on the fit seeds
%         (the model used for design),
%     (c) the -E designs predicted with the estimator-campaign pooled model that was
%         used to select them, before any of their runs was logged
%         (out of sample; both seed sets are new data for (c)).
%   Designs: PID-Budget, AA-Budget, AA-HalfGain (ids 1 2 4), PID-Budget-S
%   (new, 15), AA-Budget-S, AA-HalfGain-S (9 10), PID-E, AA-Budget-E,
%   AA-HalfGain-E (17 18 19). Stack delay convention (rates 6 ms).
%   Quantities: sigma_e (true position/yaw minus its window mean, x y z psi),
%   sigma_u and sigma_du of [u_phi u_theta u_psi u_z] (the roll/pitch
%   increments are those of the full cascade), the isolated-attitude (local)
%   roll/pitch increment prediction for reference, and J.
%   Columns: Measured_fit, Measured_heldout, prediction (a)/(b)/(c) and
%   ratios prediction/measured for each seed set.
%   Writes ../results/analysis/validation_3level.csv and noise_long/.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); DW = fullfile(px4sitl_work_dir(),'redesign_data');
D = load(fullfile(RP,'design_gains.mat')); F = load(fullfile(RA,'frozen_ceilings_and_delays.mat')); ds = F.ds;
L4 = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF0 = L4.NF;
PS = load(fullfile(RA,'pidS.mat')); E7 = load(fullfile(RA,'designs_E.mat'));
ids = [1 2 4 15 9 10 17 18 19];
des = {D.C(1), D.C(2), D.C(4), PS.PS, D.Strict(2), D.Strict(3), E7.DE(1), E7.DE(2), E7.DE(3)};
fit = 97001:97005; ho = 97006:97010; win = [20 300];
files = {};
for j = ids
    f = [dir(fullfile(DW, sprintf('LH_hover_c%d_seed*_T300.mat', j))); dir(fullfile(DW, sprintf('EFL_hover_c%d_seed*_T300.mat', j)))];
    files = [files, fullfile({f.folder}, {f.name})]; %#ok<AGROW>
end
nd = fullfile(RA,'noise_long');
if exist(fullfile(nd,'noise_models.mat'),'file'), Ln = load(fullfile(nd,'noise_models.mat')); NF = Ln.NF;
else, NF = fit_noise_models(files, ds, win, fit, ho, nd); end
names = ["se_x","se_y","se_z","se_psi","su_phi","su_theta","su_psi","su_z","sdu_phi","sdu_theta","sdu_psi","sdu_z","J"];
T = table();
for j = 1:numel(ids)
    id = ids(j); co = des{j};
    mf = measure(files, id, fit, win); mh = measure(files, id, ho, win);
    pa = vec(predict_estimator_stats(co, NF, ds, 'design', id)); pb = vec(predict_estimator_stats(co, NF, ds, 'pooled'));
    if id >= 17, pc = vec(predict_estimator_stats(co, NF0, ds, 'pooled')); else, pc = nan(1,13); end
    loc = localinc(co, NF, ds, id);
    for q = 1:13
        T = [T; table(string(co.name), id, names(q), mf.v(q), mh.v(q), pa(q), pb(q), pc(q), mf.n, mh.n, ...
            'VariableNames', {'Design','CtrlId','Quantity','Measured_fit','Measured_heldout','Pred_a_perDesign','Pred_b_pooled','Pred_c_designtime','Runs_fit','Runs_heldout'})]; %#ok<AGROW>
    end
    T = [T; table(string(co.name), id, "sdu_rp_local_pred", NaN, NaN, loc(1), loc(2), NaN, mf.n, mh.n, ...
        'VariableNames', T.Properties.VariableNames(1:10))]; %#ok<AGROW>
end
T.Ra_fit = T.Pred_a_perDesign ./ T.Measured_fit; T.Ra_ho = T.Pred_a_perDesign ./ T.Measured_heldout;
T.Rb_fit = T.Pred_b_pooled ./ T.Measured_fit;    T.Rb_ho = T.Pred_b_pooled ./ T.Measured_heldout;
T.Rc_fit = T.Pred_c_designtime ./ T.Measured_fit; T.Rc_ho = T.Pred_c_designtime ./ T.Measured_heldout;
writetable(T, fullfile(RA,'validation_3level.csv'));
end

function v = vec(z)
v = [z.se_x z.se_y z.se_z z.se_psi z.su_phi z.su_theta z.su_psi z.su_z z.sdu_phi z.sdu_theta z.sdu_psi z.sdu_z z.J];
end

function loc = localinc(co, NF, ds, id)
% isolated-attitude (local) pitch increment (class X attitude/rate errors only),
% per-design model (loc(1)) and pooled model (loc(2))
S = load(fullfile(fileparts(mfilename('fullpath')),'..','results','design_gains.mat'),'S'); S = S.S;
M = chan_model_est(co.kind, co.G(19:24)', S.Ts, [], S.act, [], [co.tauhat co.ndhat], [ds([8 11]) 0 0]);
loc = nan(1,2); L = NF.X; md = [];
for k = 1:numel(L), if numel(L(k).ctrls) == 1 && L(k).ctrls == id, md = L(k); end, end
ms = {md, NF.X_pooled};
for w = 1:2
    m = ms{w}; if isempty(m), continue, end
    m2 = m; m2.Cf = m.Cf(3:4,:); m2.Df = m.Df(3:4,:);          % attitude and rate errors only
    z = cov_estimator_aware(M, m2); loc(w) = z(3);
end
end

function m = measure(files, id, seeds, win)
% RMS over the runs of the per-run standard deviations (variances averaged)
v = []; n = 0;
for k = 1:numel(files)
    [~, nm] = fileparts(files{k}); tok = regexp(nm, '_c(\d+)_seed(\d+)_', 'tokens', 'once');
    if str2double(tok{1}) ~= id || ~ismember(str2double(tok{2}), seeds), continue, end
    D = load(files{k}); R = D.R; x = double(D.Lx); t = D.Lt;
    kt = t >= win(1) & t < win(2); se = std(x(kt, [1 2 3 9]), 1, 1);
    kw = R.tc >= win(1) & R.tc < win(2); okw = kw & R.ok; okw(1) = false; ip = find(okw);
    su = std(double(R.ua(kw,:)), 1, 1); sdu = std(double(R.ua(ip,:) - R.ua(ip-1,:)), 1, 1);
    pv = quad_params(); dW = pv(11)*double(R.wn(ip,:) - R.wn(ip-1,:)); J = sqrt(mean(dW(:).^2));
    v(end+1,:) = [se, su, sdu, J].^2; n = n + 1; %#ok<AGROW>
end
m = struct('v', sqrt(mean(v, 1)), 'n', n);
end
