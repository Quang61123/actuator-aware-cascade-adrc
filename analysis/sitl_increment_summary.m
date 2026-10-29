function [T, H, A] = sitl_increment_summary()
%SITL_INCREMENT_SUMMARY SITL summaries after the PID-Budget (b*) runs (phase P4, id 22):
%   T: PID-Budget (b*) S6/S5 RMSE, J_unq, J_q, misses, MATLAB counterparts; -E designs versus PID-E and versus
%      PID-Budget (b*): paired RMSE, J_unq, J_q (mean, CI, pairs lower); AA designs of the first campaign versus
%      PID-Budget (b*) in S6/S5; estimator-bypass (PBR) increments.
%   H: noise-only hover cascade increment per design: pooled value (larger of roll/pitch, RMS over seeds) and
%      a 95% CI over the ten seeds of the per-seed larger-axis value.
%   Writes increment_sitl.csv, increment_hover.csv.
here = fileparts(mfilename('fullpath')); RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); RR = fullfile(RA,'runs');
A = readtable(fullfile(RA,'sitl_all_runs.csv'),'TextType','string'); A.Failed = string(A.Failed)=="true" | string(A.Failed)=="1";
O = readtable(fullfile(here,'..','PX4_SITL','results','sitl_all_runs.csv'),'TextType','string'); O = O(O.TauPlant_ms==0 & O.Rep==0,:);
T = table();
add = @(T, a, b, lab, ref, qty, x, y) [T; pr(a, b, lab, ref, qty, x, y)];
for g = ["PBS","PB5","PBR"], z = A(A.Group==g,:); T = [T; table(g+" PID-Budget (b*)", "", "summary", height(z), mean(z.RMSE_mm), std(z.RMSE_mm), mean(z.J_unq), mean(z.J_q), sum(z.Miss>0), sum(z.Miss), NaN, NaN, NaN, 'VariableNames', vn())]; end %#ok<AGROW>
for gE = ["EFS","EF5"], gP = "PBS"; if gE == "EF5", gP = "PB5"; end
  P = A(A.Group==gP,:); E17 = A(A.Group==gE & A.Ctrl==17,:);
  for c = 17:19, X = A(A.Group==gE & A.Ctrl==c,:);
    for q = ["RMSE_mm","J_unq","J_q"], T = [T; pr(X, P, gE+" c"+c, "PID-Budget (b*)", q)]; if c > 17, T = [T; pr(X, E17, gE+" c"+c, "PID-E", q)]; end, end %#ok<AGROW>
  end
end
for sc = [6 5], gP = "PBS"; if sc == 5, gP = "PB5"; end; P = A(A.Group==gP,:);
  for c = ["AA-Budget","AA-HalfGain"], X = O(O.Ctrl==c & O.Scen==sc,:); X.J_unq = nan(height(X),1); X.J_q = X.Increment; T = [T; pr(X, P, "first campaign S"+sc+" "+c, "PID-Budget (b*)", "RMSE_mm")]; end %#ok<AGROW>
end
writetable(T, fullfile(RA,'increment_sitl.csv')); disp(T)
% ---- hover cascade increments ----
H = table(); src = {"PID-Budget", fullfile(RP,'estimator_noise','runs'), 'est_hover_c1_seed%d.mat'; "PID-E", RR, 'E0H_hover_c17_seed%d_T60.mat'; "AA-Budget-E", RR, 'E0H_hover_c18_seed%d_T60.mat'; "AA-HalfGain-E", RR, 'E0H_hover_c19_seed%d_T60.mat'};
for k = 1:4, V = nan(10,2);
  for q = 1:10, fn = fullfile(src{k,2}, sprintf(src{k,3}, 96000+q)); if ~exist(fn,'file'), ff = dir(fullfile(src{k,2}, strrep(sprintf(src{k,3}, 96000+q), '_T60', '*'))); if isempty(ff), continue, end, fn = fullfile(ff(1).folder, ff(1).name); end
    Z = load(fn); V(q,:) = Z.X.sd_du(1:2); end
  ok = all(isfinite(V),2); pooled = max(sqrt(mean(V(ok,:).^2,1))); ps = max(V(ok,:),[],2); n = sum(ok);
  H = [H; table(src{k,1}, n, pooled, mean(ps), tinv(.975,n-1)*std(ps)/sqrt(n), min(ps), max(ps), 0.0065137205250234, 'VariableNames', {'Design','N','Pooled','MeanPerSeed','CI','Min','Max','Bcas'})]; %#ok<AGROW>
end
H.PooledOverBcasPct = 100*(H.Pooled./H.Bcas - 1); writetable(H, fullfile(RA,'increment_hover.csv')); disp(H)
end
function v = vn(), v = {'Item','Reference','Quantity','N','Mean','SD','J_unq','J_q','RunsWithMiss','Steps','Diff','CI','Lower'}; end
function t = pr(X, P, lab, ref, q)
  [~, ia, ib] = intersect(X.Seed, P.Seed); x = X.(q)(ia); p = P.(q)(ib); d = x - p; ok = isfinite(d); d = d(ok); n = numel(d);
  t = table(lab, ref, q, n, mean(x(ok)), std(x(ok)), NaN, NaN, NaN, NaN, mean(d), tinv(.975,max(n-1,1))*std(d)/sqrt(n), sum(d<0), 'VariableNames', vn());
end