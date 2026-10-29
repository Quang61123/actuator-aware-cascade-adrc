function T = matlab_design_runs(which)
%MATLAB_DESIGN_RUNS MATLAB S6 runs of the redesigned designs.
%   which: cell array of 'A1' (A1-Commanded re-searched under the response
%   requirement, commanded_input_scan) and/or 'PIDS' (new PID-Budget-S,
%   design_pid_strict). Same protocol as mc_trials: scenario S6, 40 s,
%   window 10-40 s, seeds 91001-91030 (same exogenous realizations), true
%   rotor lag 30 ms, delay 5 samples, nominal noise. Per run: position RMSE,
%   max error, yaw RMS error, rotor-command increment, rotor fluctuation,
%   saturation, failure. Writes ../results/analysis/matlab_<which>.csv
%   and the design structs to designs_main.mat.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; P = sim_config();
fd = fullfile(RA,'designs_main.mat'); if exist(fd,'file'), M = load(fd); M = M.M; else, M = struct(); end
T = table();
for w = 1:numel(which)
    switch which{w}
        case 'A1'
            K = load(fullfile(RA,'commanded_scan.mat')); A = K.T(K.T.Design == "A1-Commanded", :);
            cells = cell(1,4); for c = 1:4, r = A(A.Class == c, :); cells{c} = [r.Best_wc r.Best_wo 0 0 0 0]; end
            co = D.C(5); co.G = pack_gains(cells); co.name = 'A1-Commanded';
            co.note = 're-searched under the response requirement'; M.A1 = co; cap = D.cap;
        case {'A1x', 'ROx'}   % ratio range extended to 1.0 (commanded_scan_ratio1)
            K = load(fullfile(RA,'commanded_scan_ratio1.mat'));
            if strcmp(which{w}, 'A1x'), nm = "A1-Commanded"; co = D.C(5); else, nm = "RO-HalfGain"; co = D.C(3); end
            A = K.T(string(K.T.Design) == nm, :);
            cells = cell(1,4); for c = 1:4, r = A(A.Class == c, :); cells{c} = [r.Best_wc r.Best_wo 0 0 0 0]; end
            co.G = pack_gains(cells); co.name = char(nm + " (ratio >= 1.0)");
            co.note = 'commanded-input design, ratio range extended to 1.0'; M.(which{w}) = co; cap = D.cap;
        case {'CURVE050', 'CURVE075'}   % budget-curve PID points re-searched (design_pid_budget_curve, OS target 9.99)
            K = load(fullfile(RA,'pid_curve_os999.mat')); gm = 0.5 + 0.25*strcmp(which{w}, 'CURVE075');
            co = K.PC(abs([K.PC.Gamma] - gm) < 1e-9).Controller; M.(which{w}) = co; cap = D.cap;
            V = verify_design(co, S, gm*D.dub, [15 5 2 5], cap); writetable(V, fullfile(RA, ['verify_' which{w} '.csv']));
        case 'PIDS'
            K = load(fullfile(RA,'pidS.mat')); co = K.PS; M.PIDS = co; cap = D.strict;
    end
    if ~startsWith(which{w}, 'CURVE'), V = verify_design(co, S, D.dub, [15 5 2 5], cap); writetable(V, fullfile(RA, ['verify_' which{w} '.csv'])); end
    Z = nan(30, 9);
    parfor q = 1:30
        seed = P.nominalSeed + q; pv = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat);
        E = make_exo(6, seed, 40, S.Ts, S); nd = 5*ones(round(40/S.Ts),1);
        if exist('run_quad_stress_mex','file') == 3, L = run_quad_stress_mex(co.kind, co.G, 6, pv, pc, 40, S.Ts, E, S.sig, nd);
        else, L = run_quad_stress(co.kind, co.G, 6, pv, pc, 40, S.Ts, E, S.sig, nd); end
        v = run_metrics(L, [10 40], pv);
        Z(q,:) = [seed, v.rmse3*1000, v.max3*1000, v.rmsepsi, v.chat_cmd, v.rotor_std, v.sat_rotor*100, v.sat_virt*100, v.diverged];
    end
    R = array2table(Z, 'VariableNames', {'Seed','RMSE_mm','MaxError_mm','YawRMS_rad','Increment','RotorStd','RotorClip_pct','VirtualClip_pct','Failed'});
    R.Design = repmat(string(co.name), 30, 1); R.Which = repmat(string(which{w}), 30, 1);
    writetable(R, fullfile(RA, ['matlab_' which{w} '.csv'])); T = [T; R]; %#ok<AGROW>
end
save(fd, 'M');
end
