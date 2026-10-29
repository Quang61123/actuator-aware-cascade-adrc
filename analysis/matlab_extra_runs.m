function matlab_extra_runs()
%MATLAB_EXTRA_RUNS MATLAB runs needed to update the summary tables for the
%   redesigned designs (same protocols as mc_trials / numerical_checks):
%   (1) PID-Budget-S (new): noise-only hover, seeds 96001-96010, 60 s, window
%       20-60 s, no bias/wind: measured roll/pitch increment and variance ratios
%       against the full-cascade prediction.
%   (2) RO-HalfGain (ratio >= 1.0, matlab_design_runs 'ROx'): S5, seeds 91001-91030;
%       noise/latency stress E5-1..6, seeds 93001-93010.
%   (3) A1-Commanded (ratio >= 1.0): CI of the S6 runs is computed from
%       matlab_A1x.csv (no new run).
%   Writes ../results/analysis/hover_pidS.csv, matlab_ROx_S5.csv,
%   matlab_ROx_stress.csv.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); P = sim_config();
D = load(fullfile(RP,'design_gains.mat')); S = D.S; pv = quad_params();
M = load(fullfile(RA,'designs_main.mat')); M = M.M; PS = load(fullfile(RA,'pidS.mat')); PS = PS.PS;
% ---- (1) hover, new PID-S ----
c = PS; Rv = nan(10, 4);
parfor q = 1:10
    E = make_exo(2, 96000+q, 60, S.Ts, S); E.bias(:) = 0; E.wind(:) = 0;
    if exist('run_quad_hover_mex','file') == 3, L = run_quad_hover_mex(c.kind, c.G, 2, pv, quad_params([], c.tauhat, c.ndhat), 60, S.Ts, E, S.sig, 5*ones(60000,1));
    else, L = run_quad_hover(c.kind, c.G, 2, pv, quad_params([], c.tauhat, c.ndhat), 60, S.Ts, E, S.sig, 5*ones(60000,1)); end
    v = L.u(L.t >= 20, :);
    Rv(q,:) = [var(v(:,4),1), var(diff(v(:,4)),1), var(v(:,5),1), var(diff(v(:,5)),1)];
end
oa = [c.tauhat c.ndhat]; Mm = chan_model_inner(c.kind, c.G(1:6)', S.Ts, c.G(19:24)', S.act, [], oa); k = Mm.keep; F = Mm.A(k,k); Bn = Mm.B(k,2:5);
cu = Mm.Ci(k); du = Mm.Di(2:5); W = diag([S.sy(1)^2 S.sv(1)^2 S.sy(3)^2 S.sv(3)^2]); Pn = dlyap(F, Bn*W*Bn');
vu = cu*Pn*cu' + du*W*du'; h = cu*(F - eye(numel(k))); qq = cu*Bn - du; vdu = h*Pn*h' + qq*W*qq' + du*W*du';
H = table(string(c.name), sqrt(vdu), sqrt(mean(Rv(:,2))), sqrt(mean(Rv(:,4))), mean(Rv(:,1))/vu, mean(Rv(:,2))/vdu, D.strict, ...
    'VariableNames', {'Ctrl','PredSigmaDU','MeasuredSigmaDU_roll','MeasuredSigmaDU_pitch','VarianceRatioU','VarianceRatioDU','CascadeCap'});
writetable(H, fullfile(RA,'hover_pidS.csv'));
% ---- (2) RO-HalfGain (ratio >= 1.0): S5 and stress ----
co = M.ROx; Z = nan(30, 9);
parfor q = 1:30
    seed = P.nominalSeed + q; pvt = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat);
    E = make_exo(5, seed, 45, S.Ts, S); nd = 5*ones(round(45/S.Ts),1);
    if exist('run_quad_s5_mex','file') == 3, L = run_quad_s5_mex(co.kind, co.G, 5, pvt, pc, 45, S.Ts, E, S.sig, nd);
    else, L = run_quad_stress(co.kind, co.G, 5, pvt, pc, 45, S.Ts, E, S.sig, nd); end
    v = run_metrics(L, [5 45], pvt);
    Z(q,:) = [seed, v.rmse3*1000, v.max3*1000, v.rmsepsi, v.chat_cmd, v.rotor_std, v.sat_rotor*100, v.sat_virt*100, v.diverged];
end
T5 = array2table(Z, 'VariableNames', {'Seed','RMSE_mm','MaxError_mm','YawRMS_rad','Increment','RotorStd','RotorClip_pct','VirtualClip_pct','Failed'});
writetable(T5, fullfile(RA,'matlab_ROx_S5.csv'));
scale = [1 .5 2 4 1 1]; pattern = {'fixed','fixed','fixed','fixed','sine','burst'}; ZS = nan(60, 6);
parfor j = 1:60
    z = ceil(j/10); q = mod(j-1,10) + 1; seed = P.stressSeed + q; pvt = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat);
    E = make_exo(6, seed, 40, S.Ts, S); n = round(40/S.Ts); t = (0:n-1)'*S.Ts; nd = 5*ones(n,1);
    if strcmp(pattern{z},'sine'), nd = round(5+5*sin(2*pi*t)); elseif strcmp(pattern{z},'burst'), nd(t >= 10 & mod(t-10,2) < .2) = 20; end
    if exist('run_quad_stress_mex','file') == 3, L = run_quad_stress_mex(co.kind, co.G, 6, pvt, pc, 40, S.Ts, E, S.sig*scale(z), nd);
    else, L = run_quad_stress(co.kind, co.G, 6, pvt, pc, 40, S.Ts, E, S.sig*scale(z), nd); end
    v = run_metrics(L, [10 40], pvt);
    ZS(j,:) = [z, seed, v.rmse3*1000, v.chat_cmd, v.rotor_std, v.diverged];
end
TS = array2table(ZS, 'VariableNames', {'Condition','Seed','RMSE_mm','Increment','RotorStd','Failed'});
TS.Label = string(arrayfun(@(z) sprintf('E5-%d', z), TS.Condition, 'uni', 0));
writetable(TS, fullfile(RA,'matlab_ROx_stress.csv'));
end
