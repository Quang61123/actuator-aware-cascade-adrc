function Rep = sitl_report_main()
%SITL_REPORT_MAIN Summary tables and S6 figure of the main PX4 SITL campaign.
%   MATLAB references: results/trials.csv (same seeds, gains,
%   plant and scenarios). Outputs: results/*.csv, results/sitl_report_main.mat,
%   ../figures/fig_sitl_s6.pdf.
here = fileparts(mfilename('fullpath'));
simRoot = fullfile(here,'..','MATLAB_simulation'); addpath(here, simRoot);
out = fullfile(here,'results'); figd = fullfile(here,'..','figures'); if ~exist(figd,'dir'), mkdir(figd); end
G = load(fullfile(here,'sitl_gain_table.mat')); names = G.names;
T = sitl_collect(); Rep.T = T;
ML = readtable(fullfile(here,'..','results','trials.csv'));
nom = T(T.TauPlant_ms == 0 & T.Rep == 0, :);
% ---------- S6 (E1, E0) and S5 / ablation summaries ----------
rows = {6,1:4,'E1'; 6,8:10,'E0'; 5,1:4,'S5'; 6,5:7,'E3'};
TS = table();
for k = 1:size(rows,1)
    s = rows{k,1}; ex = rows{k,3};
    for c = rows{k,2}
        a = sortrows(nom(nom.Scen == s & nom.id == c, :), 'Seed');
        if isempty(a), continue; end
        m = ML(strcmp(ML.Experiment, ex) & strcmp(ML.Ctrl, names{c}) & ismember(ML.Seed, a.Seed), :);
        m = sortrows(m, 'Seed');
        ok = ~a.Failed; okm = ~m.Failed;
        TS = [TS; table(string(ex), s, string(names{c}), height(a), sum(a.Failed), ...
            mean(m.RMSE_mm(okm)), std(m.RMSE_mm(okm)), mean(a.RMSE_mm(ok)), std(a.RMSE_mm(ok)), ...
            tinv(.975, sum(ok)-1)*std(a.RMSE_mm(ok))/sqrt(sum(ok)), ...
            mean(m.Increment(okm)), mean(a.Increment(ok)), mean(m.RotorStd(okm)), mean(a.RotorStd(ok)), ...
            mean(m.MaxError_mm(okm)), mean(a.Max_mm(ok)), ...
            'VariableNames',{'Exp','Scen','Ctrl','N','Fail_SITL','RMSE_ML','SD_ML','RMSE_SITL','SD_SITL','CI_SITL', ...
            'Inc_ML','Inc_SITL','RotorStd_ML','RotorStd_SITL','Max_ML','Max_SITL'})]; %#ok<AGROW>
    end
end
writetable(TS, fullfile(out,'sitl_vs_matlab_summary.csv')); Rep.TS = TS;
% ---------- paired differences (SITL and MATLAB, same seeds) ----------
pairs = [2 1 6; 4 2 6; 9 8 6; 10 8 6; 10 9 6; 2 1 5; 4 2 5];
TP = table();
for k = 1:size(pairs,1)
    a = pairs(k,1); b = pairs(k,2); s = pairs(k,3);
    A = sortrows(nom(nom.Scen == s & nom.id == a,:),'Seed'); B = sortrows(nom(nom.Scen == s & nom.id == b,:),'Seed');
    [sd, ia, ib] = intersect(A.Seed, B.Seed); if isempty(sd), continue; end
    ok = ~A.Failed(ia) & ~B.Failed(ib); d = A.RMSE_mm(ia(ok)) - B.RMSE_mm(ib(ok)); n = numel(d);
    TP = [TP; table(string(names{a}), string(names{b}), s, n, mean(d), tinv(.975,n-1)*std(d)/sqrt(n), sum(d < 0), ...
        'VariableNames',{'A','B','Scen','N','Diff_mm','CI95_mm','A_better'})]; %#ok<AGROW>
end
writetable(TP, fullfile(out,'sitl_paired.csv')); Rep.TP = TP;
% ---------- rotor-lag mismatch ----------
E = T(T.TauPlant_ms > 0, :); TE = table();
for tm = unique(E.TauPlant_ms)'
    for c = [1 2 4]
        a = E(E.TauPlant_ms == tm & E.id == c, :); if isempty(a), continue; end
        m = ML(strcmp(ML.Experiment,'E7') & strcmp(ML.Ctrl, names{c}) & abs(ML.Tau - tm/1000) < 1e-9 & ML.Delay == 5 & ismember(ML.Seed, a.Seed), :);
        TE = [TE; table(tm, string(names{c}), height(a), sum(m.Failed), sum(a.Failed), sum(m.FluctuationThreshold), sum(a.Threshold), ...
            mean(m.RMSE_mm(~m.Failed)), mean(a.RMSE_mm(~a.Failed)), ...
            'VariableNames',{'TauPlant_ms','Ctrl','N','Fail_ML','Fail_SITL','Thr_ML','Thr_SITL','RMSE_ML','RMSE_SITL'})]; %#ok<AGROW>
    end
end
writetable(TE, fullfile(out,'sitl_mismatch.csv')); Rep.TE = TE;
% ---------- repeatability, misses ----------
r = T(T.id == 2 & T.Scen == 6 & T.Seed == 91001 & T.TauPlant_ms == 0, :);
Rep.repeat = r.RMSE_mm'; Rep.missRuns = sum(T.Miss > 0); Rep.missTotal = sum(T.Miss); Rep.nRuns = height(T);
Rep.repDev = local_repdev(here);
Rep.est = local_estimation();
local_fig(here, figd, names);
save(fullfile(out,'sitl_report_main.mat'),'Rep');
end

function dv = local_repdev(here) %#ok<INUSD>
F = fullfile(px4sitl_work_dir(),'full_main'); dv = NaN;
f0 = fullfile(F,'c2_s6_seed91001_tau0_rep0.mat'); if ~exist(f0,'file'), return; end
L0 = load(f0); dv = 0;
for k = 1:2
    f = fullfile(F, sprintf('c2_s6_seed91001_tau0_rep%d.mat',k));
    if exist(f,'file'), L = load(f); n = min(size(L0.L.x,1), size(L.L.x,1));
        dv = max(dv, 1e3*max(vecnorm(L.L.x(1:n,1:3) - L0.L.x(1:n,1:3), 2, 2))); end
end
end

function local_fig(here, figd, names)
% S6, seed 91002: SITL (top) and MATLAB rerun with the same plant (bottom)
R = fullfile(here,'results','runs'); ids = [1 2 4 9];
cols = {[0.85 0.33 0.10],[0 0.45 0.74],[0.47 0.67 0.19],[0.49 0.18 0.56]};
D0 = load(fullfile(here,'..','results','design_gains.mat'));
src = [num2cell(D0.C(1:7)), num2cell(D0.Strict(1:3))]; S = D0.S;
fig = figure('Units','centimeters','Position',[2 2 8.8 8.2],'Color','w','Visible','off');
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(tl); hold(ax1,'on'); ax2 = nexttile(tl); hold(ax2,'on');
for j = 1:numel(ids)
    f = fullfile(R, sprintf('c%d_s6_seed91002_tau0_rep0.mat', ids(j))); if ~exist(f,'file'), continue; end
    Q = load(f); e = vecnorm(double(Q.D.r(:,1:3)) - double(Q.D.x(:,1:3)), 2, 2)*1e3;
    lab = strrep(names{ids(j)},'-Strict','-S');
    plot(ax1, Q.D.t, e, '-', 'Color', cols{j}, 'LineWidth', 0.9, 'DisplayName', lab);
    c = src{ids(j)}; pv = quad_params([], .03, 5); pc = quad_params([], c.tauhat, c.ndhat);
    Ex = make_exo(6, 91002, 40, S.Ts, S);
    L = run_quad_stress(c.kind, c.G, 6, pv, pc, 40, S.Ts, Ex, S.sig, 5*ones(40000,1));
    em = vecnorm(L.r(:,1:3) - L.x(:,1:3), 2, 2)*1e3;
    plot(ax2, L.t(1:10:end), em(1:10:end), '-', 'Color', cols{j}, 'LineWidth', 0.9, 'DisplayName', lab);
end
for ax = [ax1 ax2]
    xlim(ax,[0 40]); set(ax,'YScale','log','FontSize',7); ylim(ax,[1 3e4]); yticks(ax,[1 10 100 1000]);
    grid(ax,'on'); box(ax,'on'); ylabel(ax,'||e_p|| [mm]');
    legend(ax,'Location','north','NumColumns',4,'FontSize',6,'Box','off');
end
title(ax1,'(a) PX4 SITL','FontSize',7,'FontWeight','normal'); title(ax2,'(b) MATLAB simulation','FontSize',7,'FontWeight','normal');
xlabel(ax2,'t [s]');
exportgraphics(fig, fullfile(figd,'fig_sitl_s6.pdf'), 'ContentType','vector'); close(fig);
end

function E = local_estimation()
E = struct('pos_rms_mm',NaN,'att_rms_mrad',NaN,'lag_ms',NaN);
W = px4sitl_work_dir(); tag = 'c2_s6_seed91002_tau0_rep0';
uf = fullfile(W,[tag '.ulg']); lf = fullfile(W,'full_main',[tag '.mat']);
if ~exist(uf,'file') || ~exist(lf,'file'), return; end
S = load(lf); L = S.L; u = ulogreader(uf);
d = readTopicMsgs(u,'TopicNames',{'vehicle_local_position','vehicle_attitude'});
P = d.TopicMessages{strcmp(string(d.TopicNames),'vehicle_local_position')};
Q = d.TopicMessages{strcmp(string(d.TopicNames),'vehicle_attitude')};
tp = seconds(P.timestamp) - L.tArm; ok = tp >= 10 & tp < 40;
xt = interp1(L.t, L.x, tp(ok));
ep = [P.x(ok) - xt(:,1), -P.y(ok) - xt(:,2), -P.z(ok) - xt(:,3)];
E.pos_rms_mm = sqrt(mean(sum(ep.^2,2)))*1e3;
tq = seconds(Q.timestamp) - L.tArm; okq = tq >= 10 & tq < 40; q = double(Q.q(okq,:));
eul = [atan2(2*(q(:,1).*q(:,2)+q(:,3).*q(:,4)), 1-2*(q(:,2).^2+q(:,3).^2)), ...
       -asin(max(-1,min(1,2*(q(:,1).*q(:,3)-q(:,4).*q(:,2))))), ...
       -atan2(2*(q(:,1).*q(:,4)+q(:,2).*q(:,3)), 1-2*(q(:,3).^2+q(:,4).^2))];
xa = interp1(L.t, L.x(:,7:9), tq(okq)); E.att_rms_mrad = sqrt(mean(sum((eul - xa).^2,2)))*1e3;
lags = 0:30; cst = zeros(size(lags));
for k = 1:numel(lags), xs = interp1(L.t, L.x(:,1), tp(ok) - lags(k)*1e-3); cst(k) = mean((P.x(ok) - xs).^2); end
[~, k] = min(cst); E.lag_ms = lags(k);
end
