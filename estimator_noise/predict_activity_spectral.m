function T = predict_activity_spectral()
%PREDICT_ACTIVITY_SPECTRAL Spectral (non-parametric) cross-check of the estimator-aware
%   prediction. Instead of a fitted shaping filter, the
%   measured cross-spectral density matrix S(f) of the class error vector
%   (Welch, training seeds 96001-96005 of each design) is propagated through
%   the frequency response T(f) of the implemented cascade model from the
%   noise inputs to the inner virtual-control increment:
%       sigma_du^2 = integral_0^{fs/2} T(f) S(f) T(f)^H df,
%   with the same signal delays as the parametric model. Check: the same
%   integral with T = 1 - exp(-j 2 pi f Ts) on a single output must return
%   its measured increment variance (column SelfCheck).
%   Writes est_prediction_np.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RES = fullfile(RP,'estimator_noise');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; Ts = S.Ts; act = S.act;
L = load(fullfile(RES,'est_noise_models.mat')); dsig = L.NF.dsig;
DW = fullfile(px4sitl_work_dir(),'est_data');
des = {1, D.C(1); 2, D.C(2); 4, D.C(4); 8, D.Strict(1); 9, D.Strict(2); 10, D.Strict(3)};
cls = {[1 4 8 11], [1 1 1 1]; [2 5 7 10], [1 1 -1 -1]};
nf = 8192; T = table();
for r = 1:size(des,1)
    id = des{r,1}; co = des{r,2}; oa = [co.tauhat co.ndhat]; sg = zeros(1,2); chk = zeros(1,2);
    for c = 1:2
        % measured cross-spectral matrix (one-sided), averaged over segments
        Sm = 0; nn = 0; vd = 0; nd = 0;
        for s = 96001:96005
            E = estimation_errors(fullfile(DW, sprintf('est_hover_c%d_seed%d.mat', id, s)), dsig, [20 60]);
            y = E.e(:, cls{c,1}).*cls{c,2}; y = y - mean(y,1);
            [P, f] = cpsd(y, y, hann(nf), nf/2, nf, 1/Ts, 'mimo');   % P(:,a,b) = cross spectrum of a and b
            Sm = Sm + P*size(y,1); nn = nn + size(y,1);
            d = diff(y(:,4)); vd = vd + sum(d.^2); nd = nd + numel(d);
        end
        Sm = Sm/nn; df = f(2) - f(1);
        M = chan_model_est(co.kind, co.G(1:6)', Ts, co.G(19:24)', act, [], oa, dsig(cls{c,1}));
        k = M.keep; F = M.A(k,k); Bn = M.B(k,2:5); ci = M.Ci(k); di = M.Di(2:5);
        z = exp(1i*2*pi*f*Ts); acc = 0; acc1 = 0;
        for q = 1:numel(f)
            G = ci/(z(q)*eye(numel(k)) - F)*Bn + di;      % noise -> inner virtual control
            Tq = (1 - 1/z(q))*G;                           % increment
            Sq = reshape(Sm(q,:,:), 4, 4);
            acc = acc + real(conj(Tq)*Sq*Tq.');
            acc1 = acc1 + abs(1 - 1/z(q))^2*real(Sq(4,4));
        end
        sg(c) = sqrt(acc*df); chk(c) = (acc1*df)/(vd/nd);
    end
    T = [T; table(string(co.name), max(sg), sg(1), sg(2), chk(1), chk(2), ...
        'VariableNames', {'Design','NP_SigmaDU','NP_SigmaDU_pitch','NP_SigmaDU_roll','SelfCheck_X','SelfCheck_Y'})]; %#ok<AGROW>
end
writetable(T, fullfile(RES,'est_prediction_np.csv'));
end
