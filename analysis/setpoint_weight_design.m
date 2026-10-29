function [cs, B] = setpoint_weight_design(co, isE)
%SETPOINT_WEIGHT_DESIGN b* rule of setpoint_rule.txt for one PID design (same code path as setpoint_weight_select):
%   per class, largest b in ({0,0.05,...,1} U {b_cur}) n [b_cur,1] with OS <= 10 % and ts <= [15 5 2 5] s;
%   order phth, z, psi, then xy on the cascade with the new phth; revert phth if xy has no feasible b.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); D = load(fullfile(RP,'design_gains.mat')); S = D.S;
F = load(fullfile(RP,'analysis','frozen_ceilings_and_delays.mat')); ds = F.ds; oa = [co.tauhat co.ndhat];
req = [15 5 2 5]; grid = 0:0.05:1; ib = [5 17 23 35]; cs = co; B = table();
for h = [3 2 4 1]
  g0 = cs.G(ib(h)); cand = unique([grid(grid >= g0) g0]); best = g0; osb = NaN; tsb = NaN;
  for b = sort(cand,'descend'), c2 = setb(cs, h, b); r = resp(c2, h, isE, oa, ds, S); if r(1) <= S.osmax + 1e-9 && r(2) <= req(h), best = b; osb = r(1); tsb = r(2); break, end, end
  if h == 1 && isnan(osb), cs = setb(cs, 3, co.G(ib(3))); for b = sort(cand,'descend'), c2 = setb(cs, 1, b); r = resp(c2, 1, isE, oa, ds, S); if r(1) <= S.osmax + 1e-9 && r(2) <= req(1), best = b; osb = r(1); tsb = r(2); break, end, end, end
  r0 = resp(co, h, isE, oa, ds, S); cs = setb(cs, h, best);
  B = [B; table(string(co.name), h, g0, best, r0(1), r0(2), osb, tsb, 'VariableNames', {'Design','Class','b_cur','b_star','OS_cur','ts_cur','OS_star','ts_star'})]; %#ok<AGROW>
end
end
function c = setb(c, h, b), ib = [5 17 23 35]; c.G(ib(h)) = b; if h == 1 || h == 3, c.G(ib(h)+6) = b; end, end
function r = resp(co, h, isE, oa, ds, S)
  ix = {1:6, 13:18, 19:24, 31:36}; g = co.G(ix{h})'; gin = []; if h == 1, gin = co.G(19:24)'; end
  if isE, r = est_stack_tools('resp', co.kind, g, h, gin, oa, ds); else, [os, ts] = step_response_metrics(chan_model(co.kind, g, S.Ts, gin, S.act, [], oa), S.Ts, 1, 400); r = [os ts]; end
end