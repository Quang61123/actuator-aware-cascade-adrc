function [T, Tchk, Tms] = predict_activity()
%PREDICT_ACTIVITY Estimator-aware activity prediction vs. PX4 SITL.
%   White-noise consistency check: with a white "filter" (Af empty, Df = noise standard deviations)
%   and the common delay (5 samples), cov_estimator_aware with
%   chan_model_est reproduces the white-noise covariances
%   (cascade_stats / chan_var) -> est_cov_check.csv.
%   Predictions: for the six designs (gains unchanged) the cascade-level
%   inner increment sigma_du (roll/pitch, larger of the two, as cascade_stats)
%   and the rotor-command increment J (linearized mixer at
%   hover, channels assumed uncorrelated) are predicted
%     (a) with the white measurement-noise model and
%     (b) with the estimator-coloured noise models of fit_noise_per_design and the
%         signal delays of estimate_delays,
%   and compared with the noise-only hover SITL runs with EKF2 input
%   (est_hover_*, seeds 96001-96010; also on the held-out seeds 96006-96010
%   that were not used for fitting) -> est_prediction.csv.
%   Rotor-formula check: the white prediction of J is compared with the
%   estimator-bypass hover runs of the activity check (act_raw_hover_*).
%   Input-sensitivity peak Ms of the three main designs with the
%   estimator delays -> est_ms.csv.
here = fileparts(mfilename('fullpath'));
simRoot = fullfile(here,'..','MATLAB_simulation'); addpath(here, simRoot);
RP = fullfile(here,'..','results'); RES = fullfile(RP,'estimator_noise');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; pv = quad_params();
L = load(fullfile(RES,'est_noise_models.mat')); NF = L.NF; dsig = NF.dsig;
des = {1, D.C(1); 2, D.C(2); 4, D.C(4); 8, D.Strict(1); 9, D.Strict(2); 10, D.Strict(3)};
Ts = S.Ts; act = S.act;
% linearized mixer at hover: dw = Jm * [dT; dtau]
k = pv(4); l = pv(3); b = pv(5); wh = sqrt(pv(1)*pv(2)/(4*k));
Mx = [k k k k; 0 -l*k 0 l*k; -l*k 0 l*k 0; -b b -b b]; Jm = diag(1/(2*wh)*ones(4,1))*inv(Mx); %#ok<MINV>
I3 = pv(7:9);
% ---------------- white-noise consistency check ----------------
Tchk = table();
for r = 1:size(des,1)
    co = des{r,2}; oa = [co.tauhat co.ndhat];
    z0 = cascade_stats(co, S);
    M0 = chan_model_inner(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa);
    M1 = chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, act(2)*[1 1 1 1]);
    same = isequal(M0.A, M1.A) && isequal(M0.B, M1.B) && isequal(M0.Ci, M1.Ci) && isequal(M0.Di, M1.Di);
    W = white(S.sy(1), S.sv(1), S.sy(3), S.sv(3)); z1 = cov_estimator_aware(M1, W);
    e1 = max(abs(z1(2:4)./z0(2:4) - 1));
    es = zeros(1,2); hh = [2 4]; gi = {13:18, 31:36};
    for q = 1:2
        g = co.G(gi{q})'; [~,~,vu,vdu] = chan_var(co.kind, g, Ts, S.sy(hh(q)), S.sv(hh(q)), 0, [], 0, 0, act, oa);
        Mq = chan_model_est(co.kind, g, Ts, [], act, [], oa, act(2)*[1 1 1 1]);
        zq = cov_estimator_aware(Mq, white(S.sy(hh(q)), S.sv(hh(q))));
        es(q) = max(abs([zq(2)^2/vu, zq(3)^2/vdu] - 1));
    end
    Tchk = [Tchk; table(string(co.name), same, e1, es(1), es(2), ...
        'VariableNames', {'Design','ModelIdentical','CascadeRelErr','ZRelErr','YawRelErr'})]; %#ok<AGROW>
end
writetable(Tchk, fullfile(RES,'est_cov_check.csv'));
% ---------------- measured SITL ----------------
RR = fullfile(RES,'runs'); RACT = fullfile(RP,'sitl_activity','runs');
% ---------------- predictions ----------------
T = table();
for r = 1:size(des,1)
    id = des{r,1}; co = des{r,2}; oa = [co.tauhat co.ndhat];
    % (a) white model (delay 5 samples)
    zw = cascade_stats(co, S);
    sw = chan_sig(co, Ts, act, oa, act(2)*[1 1 1 1], act(2)*[1 1 1 1], act(2)*[1 1], act(2)*[1 1], ...
        white(S.sy(1),S.sv(1),S.sy(3),S.sv(3)), white(S.sy(1),S.sv(1),S.sy(3),S.sv(3)), white(S.sy(2),S.sv(2)), white(S.sy(4),S.sv(4)));
    Jw = rotorJ(sw, Jm, pv(1), I3);
    % (b) estimator-coloured model
    se = chan_sig(co, Ts, act, oa, dsig([1 4 8 11]), dsig([2 5 7 10]), dsig([3 6]), dsig([9 12]), ...
        pick(NF.X, id), pick(NF.Y, id), pick(NF.Z, id), pick(NF.P, id));
    Je = rotorJ(se, Jm, pv(1), I3);
    % (b') same with the noise models pooled over all designs
    sp = chan_sig(co, Ts, act, oa, dsig([1 4 8 11]), dsig([2 5 7 10]), dsig([3 6]), dsig([9 12]), ...
        pooled(NF,'X'), pooled(NF,'Y'), pooled(NF,'Z'), pooled(NF,'P'));
    Jp = rotorJ(sp, Jm, pv(1), I3);
    % measured (EKF2 input, hover)
    [mall, Jall, nall] = measured(RR, sprintf('est_hover_c%d_seed*.mat', id), 96001:96010);
    [mval, Jval, nval] = measured(RR, sprintf('est_hover_c%d_seed*.mat', id), 96006:96010);
    [~, Jraw] = measured(RACT, sprintf('act_raw_hover_c%d_seed*.mat', id), 96001:96010);
    T = [T; table(string(co.name), zw(4), max(se.dui), mall, mall/max(se.dui), mval, mval/max(se.dui), ...
        Jw, Jraw, Je, Jall, Jall/Je, Jval, Jval/Je, nall, nval, se.dui(1), se.dui(2), se.rho, ...
        max(sp.dui), mall/max(sp.dui), Jp, Jall/Jp, ...
        'VariableNames', {'Design','Pred_white_SigmaDU','Pred_est_SigmaDU','SITL_SigmaDU','Ratio_SigmaDU', ...
        'SITL_SigmaDU_heldout','Ratio_SigmaDU_heldout','Pred_white_J','SITLraw_J','Pred_est_J','SITL_J','Ratio_J', ...
        'SITL_J_heldout','Ratio_J_heldout','N','N_heldout','Pred_est_SigmaDU_pitch','Pred_est_SigmaDU_roll','Rho', ...
        'Pred_pooled_SigmaDU','Ratio_pooled_SigmaDU','Pred_pooled_J','Ratio_pooled_J'})]; %#ok<AGROW>
end
writetable(T, fullfile(RES,'est_prediction.csv'));
% ---------------- Ms with the estimator delays ----------------
Tms = table();
for r = 1:3
    co = des{r,2}; oa = [co.tauhat co.ndhat];
    cl = {1:6, 19:24, 13:18, 19:24, 31:36};
    m0 = zeros(1,4); m1 = zeros(1,4);
    dW = {act(2)*[1 1 1 1], act(2)*[1 1 1 1], act(2)*[1 1 1 1], act(2)*[1 1 1 1]};
    dE = {dsig([1 4 8 11]), [dsig([3 6]) 0 0], [dsig([8 11]) 0 0], [dsig([9 12]) 0 0]};
    gs = {co.G(1:6)', co.G(13:18)', co.G(19:24)', co.G(31:36)'}; gin = {co.G(19:24)', [], [], []};
    for h = 1:4
        m0(h) = msens(chan_model_est(co.kind, gs{h}, Ts, gin{h}, act, [], oa, dW{h}), Ts);
        m1(h) = msens(chan_model_est(co.kind, gs{h}, Ts, gin{h}, act, [], oa, dE{h}), Ts);
    end
    Tms = [Tms; table(string(co.name), max(m0), max(m1), m0, m1, 'VariableNames', {'Design','Ms_white','Ms_est','Ms_white_classes','Ms_est_classes'})]; %#ok<AGROW>
    clear cl
end
writetable(splitvars(Tms), fullfile(RES,'est_ms.csv'));
end

% ------------------------------------------------------------------------
function W = white(varargin)
s = [varargin{:}]; p = numel(s);
W = struct('Af', zeros(0,0), 'Bf', zeros(0,p), 'Cf', zeros(p,0), 'Df', diag(s));
end

function m = pooled(NF, c)
if isfield(NF, [c '_pooled']), m = NF.([c '_pooled']); else, m = NF.(c)(1); end
end

function m = pick(list, id)
m = list(1);
for i = 1:numel(list), if any(list(i).ctrls == id), m = list(i); end, end
end

function s = chan_sig(co, Ts, act, oa, dX, dY, dZ, dP, nX, nY, nZ, nP)
% predicted increment standard deviations of all channels
zX = cov_estimator_aware(chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, dX), nX);
zY = cov_estimator_aware(chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, dY), nY);
zZ = cov_estimator_aware(chan_model_est(co.kind, co.G(13:18)', Ts, [], act, [], oa, [dZ 0 0]), nZ);
zP = cov_estimator_aware(chan_model_est(co.kind, co.G(31:36)', Ts, [], act, [], oa, [dP 0 0]), nP);
s.dui = [zX(4), zY(4)];          % inner (pitch from x, roll from y) virtual-control increments
s.duz = zZ(3); s.dupsi = zP(3); s.rho = max([zX(5) zY(5) zZ(5) zP(5)]);
end

function J = rotorJ(s, Jm, m, I3)
% J = sqrt(mean_i var(dw_i)), inputs [T tau_phi tau_theta tau_psi] uncorrelated
v = [(m*s.duz)^2, (I3(1)*s.dui(2))^2, (I3(2)*s.dui(1))^2, (I3(3)*s.dupsi)^2];
J = sqrt(mean((Jm.^2)*v(:)));
end

function [sdu, J, n] = measured(folder, pat, seeds)
f = dir(fullfile(folder, pat)); v = []; jj = [];
for i = 1:numel(f)
    R = load(fullfile(folder, f(i).name)); if ~ismember(R.run.seed, seeds), continue; end
    v(end+1,:) = R.X.sd_du(1:2).^2; jj(end+1) = R.X.J_unq^2; %#ok<AGROW>
end
n = size(v,1);
if n == 0, sdu = NaN; J = NaN; return, end
sdu = max(sqrt(mean(v,1))); J = sqrt(mean(jj));
end

function Ms = msens(M, Ts)
% input-sensitivity peak at the actuator input (as chan_ms)
k = M.keep; A = M.A(k,k); B = M.B(k,7); C = M.Ct(k); Dd = M.Dt(7);
if max(abs(eig(A))) >= 1, Ms = Inf; return, end
sys = ss(A,B,C,Dd,Ts); wg = [0 logspace(-2,log10(pi/Ts),4000)];
Ms = max(abs(squeeze(freqresp(sys,wg))));
end
