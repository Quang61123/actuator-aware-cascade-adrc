function [B, Q] = setpoint_weight_select()
%SETPOINT_WEIGHT_SELECT Largest admissible set-point weight b* per class for the PID designs,
%   following the rule declared in setpoint_rule.txt (written before this script ran). Also checks that
%   Ms, sigma_u, sigma_du, the cascade increment and (PID-E) the estimator-aware increments do not
%   depend on b. Writes setpoint_weight_select.csv, setpoint_invariance.csv and setpoint_designs.mat (struct C: cur/star).
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; K = load(fullfile(RA,'pidS.mat')); E = load(fullfile(RA,'designs_E.mat'));
L = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF = L.NF; F = load(fullfile(RA,'frozen_ceilings_and_delays.mat')); ds = F.ds;
des = {D.C(1), false; K.PS, false; E.DE(1), true}; req = [15 5 2 5]; grid = 0:0.05:1; ib = [5 13+4 19+4 31+4];  % b index per class (xy,z,phth,psi)
B = table(); Q = table(); C = struct('cur',{},'star',{});
for d = 1:3
  co = des{d,1}; isE = des{d,2}; oa = [co.tauhat co.ndhat]; cs = co;
  for h = [3 2 4 1]
    g0 = cs.G(ib(h)); cand = unique([grid(grid >= g0) g0]); best = g0; osb = NaN; tsb = NaN;
    for b = sort(cand,'descend')
      c2 = cs; c2.G(ib(h)) = b; if h == 3, c2.G(ib(h)+6) = b; end; if h == 1, c2.G(ib(h)+6) = b; end
      r = resp(c2, h, isE, oa, ds, S); if r(1) <= S.osmax + 1e-9 && r(2) <= req(h), best = b; osb = r(1); tsb = r(2); break, end
    end
    if h == 1 && isnan(osb)   % rule: revert phth to b_cur and re-evaluate xy
      cs.G([ib(3) ib(3)+6]) = co.G([ib(3) ib(3)+6]);
      for b = sort(cand,'descend'), c2 = cs; c2.G([ib(1) ib(1)+6]) = b; r = resp(c2, 1, isE, oa, ds, S); if r(1) <= S.osmax + 1e-9 && r(2) <= req(1), best = b; osb = r(1); tsb = r(2); break, end, end
    end
    r0 = resp(co, h, isE, oa, ds, S); cs.G(ib(h)) = best; if h == 1 || h == 3, cs.G(ib(h)+6) = best; end
    B = [B; table(string(co.name), h, g0, best, r0(1), r0(2), osb, tsb, 'VariableNames', {'Design','Class','b_cur','b_star','OS_cur','ts_cur','OS_star','ts_star'})]; %#ok<AGROW>
  end
  cs.name = [co.name ' (b*)']; C(d).cur = co; C(d).star = cs;
  % invariance check
  for h = 1:4
    gs0 = cls(co, h); gs1 = cls(cs, h); gin0 = []; gin1 = []; if h == 1, gin0 = co.G(19:24)'; gin1 = cs.G(19:24)'; end
    q0 = channel_stats(co.kind, gs0, h, S, gin0, oa, S.act, false, 4000); q1 = channel_stats(cs.kind, gs1, h, S, gin1, oa, S.act, false, 4000);
    v0 = q0(2:4); v1 = q1(2:4); nm = {'sigma_u','sigma_du','Ms'};
    if isE, m0 = est_stack_tools('ms', co.kind, gs0, h, gin0, oa, ds, 4000); m1 = est_stack_tools('ms', cs.kind, gs1, h, gin1, oa, ds, 4000); v0 = [v0 m0]; v1 = [v1 m1]; nm = [nm {'Ms_E'}]; end
    for k = 1:numel(v0), Q = [Q; table(string(co.name), h, string(nm{k}), v0(k), v1(k), abs(v1(k)-v0(k))/max(abs(v0(k)),eps), 'VariableNames', {'Design','Class','Quantity','Cur','Star','RelDiff'})]; end %#ok<AGROW>
  end
  if isE, z0 = est_stack_tools('cascade', co, S, NF, ds); z1 = est_stack_tools('cascade', cs, S, NF, ds); y0 = est_stack_tools('single', co.kind, cls(co,4), 4, S, NF, ds, oa); y1 = est_stack_tools('single', cs.kind, cls(cs,4), 4, S, NF, ds, oa);
    Q = [Q; table([string(co.name); string(co.name)], [1;4], ["cascade_du_est";"psi_du_est"], [z0(4); y0(3)], [z1(4); y1(3)], abs([z1(4)-z0(4); y1(3)-y0(3)])./[z0(4); y0(3)], 'VariableNames', Q.Properties.VariableNames)];
  else, z0 = cascade_stats(co, S); z1 = cascade_stats(cs, S); Q = [Q; table(string(co.name), 1, "cascade_du", z0(4), z1(4), abs(z1(4)-z0(4))/z0(4), 'VariableNames', Q.Properties.VariableNames)]; end
end
writetable(B, fullfile(RA,'setpoint_weight_select.csv')); writetable(Q, fullfile(RA,'setpoint_invariance.csv')); save(fullfile(RA,'setpoint_designs.mat'), 'C');
disp(B); fprintf('max rel diff (invariance): %.3g\n', max(Q.RelDiff));
end
function g = cls(co, h), ix = {1:6, 13:18, 19:24, 31:36}; g = co.G(ix{h})'; end
function r = resp(co, h, isE, oa, ds, S)
  g = cls(co, h); gin = []; if h == 1, gin = co.G(19:24)'; end
  if isE, r = est_stack_tools('resp', co.kind, g, h, gin, oa, ds); else, [os, ts] = step_response_metrics(chan_model(co.kind, g, S.Ts, gin, S.act, [], oa), S.Ts, 1, 400); r = [os ts]; end
end