function [R, T] = setpoint_weight_runs()
%SETPOINT_WEIGHT_RUNS MATLAB S6 (91001-91030, 40 s, window 10-40 s) and S5 (91001-91010, 45 s,
%   window 5-45 s) runs of the PID designs with the current and the b* set-point weights
%   (setpoint_designs.mat), nominal plant (tau 30 ms, delay 5). All runs reported.
%   Writes setpoint_weight_runs.csv (per run) and setpoint_summary.csv (mean, SD, paired change b* - b_cur, CI, pairs better).
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); D = load(fullfile(RP,'design_gains.mat')); S = D.S; P = sim_config();
K = load(fullfile(RA,'setpoint_designs.mat')); C = K.C; R = table();
cases = {6, 30, 40, [10 40]; 5, 10, 45, [5 45]};
for d = 1:numel(C), for v = 1:2, co = C(d).cur; if v == 2, co = C(d).star; end
  for s = 1:2, sc = cases{s,1}; n = cases{s,2}; Te = cases{s,3}; win = cases{s,4}; Z = nan(n, 6);
    parfor q = 1:n
      seed = P.nominalSeed + q; pv = quad_params([], .03, 5); pc = quad_params([], co.tauhat, co.ndhat);
      E = make_exo(sc, seed, Te, S.Ts, S); nd = 5*ones(round(Te/S.Ts),1);
      if sc == 6 && exist('run_quad_stress_mex','file') == 3, L = run_quad_stress_mex(co.kind, co.G, sc, pv, pc, Te, S.Ts, E, S.sig, nd); elseif sc == 5 && exist('run_quad_s5_mex','file') == 3, L = run_quad_s5_mex(co.kind, co.G, sc, pv, pc, Te, S.Ts, E, S.sig, nd); else, L = run_quad_stress(co.kind, co.G, sc, pv, pc, Te, S.Ts, E, S.sig, nd); end; m = run_metrics(L, win, pv);
      Z(q,:) = [seed, m.rmse3*1000, m.max3*1000, m.rmsepsi, m.chat_cmd, m.diverged];
    end
    t = array2table(Z, 'VariableNames', {'Seed','RMSE_mm','MaxError_mm','YawRMS_rad','Increment','Failed'});
    t.Design = repmat(string(C(d).cur.name), n, 1); t.Config = repmat(string(v), n, 1); t.Scen = repmat(sc, n, 1); R = [R; t]; %#ok<AGROW>
  end
end, end
R.Config(R.Config=="1") = "b_cur"; R.Config(R.Config=="2") = "b_star"; writetable(R, fullfile(RA,'setpoint_weight_runs.csv'));
T = table();
for d = 1:numel(C), for sc = [6 5]
  a = R(R.Design==string(C(d).cur.name) & R.Scen==sc & R.Config=="b_cur",:); b = R(R.Design==string(C(d).cur.name) & R.Scen==sc & R.Config=="b_star",:);
  [~,ia,ib] = intersect(a.Seed, b.Seed); x = b.RMSE_mm(ib) - a.RMSE_mm(ia); n = numel(x); tt = tinv(0.975, n-1);
  T = [T; table(string(C(d).cur.name), sc, n, mean(a.RMSE_mm), std(a.RMSE_mm), mean(b.RMSE_mm), std(b.RMSE_mm), mean(x), tt*std(x)/sqrt(n), sum(x<0), 100*mean(x)/mean(a.RMSE_mm), mean(a.Increment), mean(b.Increment), sum(a.Failed)+sum(b.Failed), ...
    'VariableNames', {'Design','Scen','N','RMSE_cur','SD_cur','RMSE_star','SD_star','Diff','CI','PairsBetter','ChangePct','Inc_cur','Inc_star','Failed'})]; %#ok<AGROW>
end, end
writetable(T, fullfile(RA,'setpoint_summary.csv')); disp(T)
end