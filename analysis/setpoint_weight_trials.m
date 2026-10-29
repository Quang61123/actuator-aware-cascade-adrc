function T = setpoint_weight_trials()
%SETPOINT_WEIGHT_TRIALS MATLAB reruns of every PID-Budget experiment whose
%   reference changes, with the b* set-point weights (setpoint_designs.mat, rule setpoint_rule.txt), using the
%   protocol of mc_trials: E1 (S6, 91001-91030), E5-1..6 (stress, 93001-93010), E7 (mismatch grid,
%   92001-92010), S5 (91001-91030); budget-curve PID points (E4): b* rule applied to every PID curve design
%   (design_gains Curve gamma >= 1, pid_curve_os999 gamma 0.5/0.75), run with b_cur and b*; noise-only hover
%   (96001-96010, 60 s, window 20-60 s) with b_cur and b* to check that hover results do not change.
%   Writes setpoint_trials.csv (per run), setpoint_curve_bstar.csv, setpoint_hover_check.csv.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); P = sim_config(); D = load(fullfile(RP,'design_gains.mat')); S = D.S;
K = load(fullfile(RA,'setpoint_designs.mat')); PB = K.C(1).star; PB.name = 'PID-Budget (b*)';
% --- curve PID designs and their b* ---
CV = struct('co',{},'gamma',{}); for j = 1:numel(D.Curve), if strcmp(D.Curve(j).Family,'PID') && D.Curve(j).Gamma >= 1, CV(end+1) = struct('co', D.Curve(j).Controller, 'gamma', D.Curve(j).Gamma); end, end %#ok<AGROW>
Q = load(fullfile(RA,'pid_curve_os999.mat')); for j = 1:numel(Q.PC), CV(end+1) = struct('co', Q.PC(j).Controller, 'gamma', Q.PC(j).Gamma); end %#ok<AGROW>
BC = table(); jobs = struct('Exp',{},'co',{},'seed',{},'tau',{},'nd',{},'noise',{},'pat',{},'scen',{},'gamma',{},'cfg',{});
for j = 1:numel(CV)
  [cs, b] = setpoint_weight_design(CV(j).co, false); b.Gamma = repmat(CV(j).gamma, height(b), 1); BC = [BC; b]; %#ok<AGROW>
  for q = 1:30, jobs(end+1) = J('E4', CV(j).co, P.nominalSeed+q, .03, 5, 1, 'fixed', 6, CV(j).gamma, 'b_cur'); jobs(end+1) = J('E4', cs, P.nominalSeed+q, .03, 5, 1, 'fixed', 6, CV(j).gamma, 'b_star'); end %#ok<AGROW>
end
writetable(BC, fullfile(RA,'setpoint_curve_bstar.csv'));
for q = 1:30, jobs(end+1) = J('E1', PB, P.nominalSeed+q, .03, 5, 1, 'fixed', 6, NaN, 'b_star'); jobs(end+1) = J('S5', PB, P.nominalSeed+q, .03, 5, 1, 'fixed', 5, NaN, 'b_star'); end %#ok<AGROW>
scale = [1 .5 2 4 1 1]; pat = {'fixed','fixed','fixed','fixed','sine','burst'};
for z = 1:6, for q = 1:10, jobs(end+1) = J(sprintf('E5-%d',z), PB, P.stressSeed+q, .03, 5, scale(z), pat{z}, 6, NaN, 'b_star'); end, end %#ok<AGROW>
for tau = P.taus, for nd = P.delays, for q = 1:10, jobs(end+1) = J('E7', PB, P.gridSeed+q, tau, nd, 1, 'fixed', 6, NaN, 'b_star'); end, end, end %#ok<AGROW>
n = numel(jobs); Z = nan(n, 7);
parfor k = 1:n
  j = jobs(k); pv = quad_params([], j.tau, j.nd); pc = quad_params([], j.co.tauhat, j.co.ndhat); dur = 40; win = [10 40]; if j.scen == 5, dur = 45; win = [5 45]; end
  E = make_exo(j.scen, j.seed, dur, S.Ts, S); nn = round(dur/S.Ts); t = (0:nn-1)'*S.Ts; ndv = j.nd*ones(nn,1);
  if strcmp(j.pat,'sine'), ndv = round(5+5*sin(2*pi*t)); elseif strcmp(j.pat,'burst'), ndv(t>=10 & mod(t-10,2)<.2) = 20; end
  if j.scen ~= 5 && exist('run_quad_stress_mex','file') == 3, L = run_quad_stress_mex(j.co.kind, j.co.G, j.scen, pv, pc, dur, S.Ts, E, S.sig*j.noise, ndv);
  elseif j.scen == 5 && exist('run_quad_s5_mex','file') == 3, L = run_quad_s5_mex(j.co.kind, j.co.G, j.scen, pv, pc, dur, S.Ts, E, S.sig*j.noise, ndv);
  else, L = run_quad_stress(j.co.kind, j.co.G, j.scen, pv, pc, dur, S.Ts, E, S.sig*j.noise, ndv); end
  v = run_metrics(L, win, pv); Z(k,:) = [v.rmse3*1000, v.chat_cmd, v.rotor_std, v.sat_rotor*100, v.sat_virt*100, v.diverged, v.max3*1000];
end
T = table(string({jobs.Exp})', string(cellfun(@(c) c.name, {jobs.co}, 'uni', 0))', [jobs.seed]', [jobs.tau]', [jobs.nd]', [jobs.noise]', string({jobs.pat})', [jobs.scen]', [jobs.gamma]', string({jobs.cfg})', ...
  'VariableNames', {'Experiment','Ctrl','Seed','Tau','Delay','Noise','Pattern','Scenario','Gamma','Config'});
T.RMSE_mm = Z(:,1); T.Increment = Z(:,2); T.RotorStd = Z(:,3); T.RotorClip_pct = Z(:,4); T.VirtualClip_pct = Z(:,5); T.Failed = Z(:,6) > 0; T.MaxError_mm = Z(:,7); T.FluctuationThreshold = Z(:,3) > 20 & ~T.Failed;
writetable(T, fullfile(RA,'setpoint_trials.csv'));
% --- hover check (b_cur vs b*) ---
pv = quad_params(); H = nan(10, 9); C0 = K.C(1).cur; C1 = K.C(1).star;
parfor q = 1:10
  E = make_exo(2, 96000+q, 60, S.Ts, S); E.bias(:) = 0; E.wind(:) = 0; out = zeros(2,4);
  for v = 1:2, c = C0; if v == 2, c = C1; end
    if exist('run_quad_hover_mex','file') == 3, L = run_quad_hover_mex(c.kind, c.G, 2, pv, quad_params([], c.tauhat, c.ndhat), 60, S.Ts, E, S.sig, 5*ones(60000,1));
    else, L = run_quad_hover(c.kind, c.G, 2, pv, quad_params([], c.tauhat, c.ndhat), 60, S.Ts, E, S.sig, 5*ones(60000,1)); end
    u = L.u(L.t >= 20, :); out(v,:) = [std(diff(u(:,4)),1), std(diff(u(:,5)),1), std(u(:,4),1), std(u(:,5),1)];
  end
  H(q,:) = [96000+q, out(1,:), out(2,:)];
end
H = array2table(H, 'VariableNames', {'Seed','sdDu_roll_cur','sdDu_pitch_cur','sdU_roll_cur','sdU_pitch_cur','sdDu_roll_star','sdDu_pitch_star','sdU_roll_star','sdU_pitch_star'});
H.MaxRelDiff = max(abs([H{:,6:9}] - [H{:,2:5}])./abs([H{:,2:5}]), [], 2); writetable(H, fullfile(RA,'setpoint_hover_check.csv'));
fprintf('hover max rel diff %.3g\n', max(H.MaxRelDiff));
end
function j = J(e, co, seed, tau, nd, noise, pat, scen, gamma, cfg)
j = struct('Exp', e, 'co', co, 'seed', seed, 'tau', tau, 'nd', nd, 'noise', noise, 'pat', pat, 'scen', scen, 'gamma', gamma, 'cfg', cfg);
end