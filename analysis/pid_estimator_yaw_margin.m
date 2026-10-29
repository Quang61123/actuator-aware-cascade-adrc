function pid_estimator_yaw_margin()
%PID_ESTIMATOR_YAW_MARGIN PID-E yaw search repeated with an overshoot design target of
%   9.99 % (the SQP solution of design_estimator_stack meets OS <= 10 % only to the solver
%   tolerance, and neither the set-point weight nor bandwidth scaling brings it
%   below 10 %). All other constraints and the 20 starts are unchanged; the
%   final check is against the specification (OS <= 10 %, 4000 frequencies).
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub;
L = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF = L.NF; F = load(fullfile(RA,'frozen_ceilings_and_delays.mat'));
g0 = D.C(1).G(31:36)'; oa = [D.C(1).tauhat D.C(1).ndhat]; Sd = S; Sd.osmax = 9.99; rng(70003,'twister'); n = 20;
starts = [{zeros(1,5)}, arrayfun(@(k) 0.6*randn(1,5), 2:n, 'uni', 0)]; res = cell(1,n);
parfor k = 1:n, res{k} = fit_pid_yaw(g0, starts{k}, Sd, dub, NF, F.ds, oa, F.Bpsi, 900); end
T = table(); best = []; bj = inf;
for k = 1:n
 r = res{k}; rs = est_stack_tools('resp', 2, r.g, 4, [], oa, F.ds); ok = r.feasible && rs(1) <= S.osmax;
 T = [T; table(k, r.J, r.feasible, rs(1), rs(2), 'VariableNames', {'Start','sigma_e','FeasibleAtTarget','OS','ts'})]; %#ok<AGROW>
 if ok && r.J < bj, bj = r.J; best = r.g; end
end
writetable(T, fullfile(RA,'pidE_yaw_margin_starts.csv')); save(fullfile(RA,'pid_estimator_yaw_margin.mat'), 'best', 'T');
end
