function [Tsitl, Thover, Ttrack] = regenerate_tables(root)
%REGENERATE_TABLES Regenerates the SITL/stack summary tables (sitl_vs_matlab, stack_hover_table,
%   stack_tracking_table) from CSV/MAT result files
%   only (no simulation), so that they can be checked from a clean copy of the results folder.
%   root: folder that contains results/ and PX4_SITL/results/ (default: parent of this file's folder).
%   Inputs : PX4_SITL/results/sitl_vs_matlab_summary.csv, results/analysis/{setpoint_trials, setpoint_weight_runs,
%            matlab_PIDS, matlab_ROx, matlab_ROx_S5, matlab_A1x, sitl_all_runs,
%            designs_E_verify, check_nominal_on_stack, est_design_iteration0, sitl_hover_E, increment_hover, increment_sitl}.csv
%   Outputs: results/analysis/sitl_vs_matlab.csv, stack_hover_table.csv, stack_tracking_table.csv
if nargin < 1, root = fullfile(fileparts(mfilename('fullpath')), '..'); end
RA = fullfile(root,'results','analysis'); rd = @(f) readtable(fullfile(RA,f),'TextType','string');
V = readtable(fullfile(root,'PX4_SITL','results','sitl_vs_matlab_summary.csv'),'TextType','string');
A = rd('sitl_all_runs.csv'); C1 = rd('setpoint_trials.csv'); C1r = rd('setpoint_weight_runs.csv');
ms = @(x) [mean(x) std(x)]; first10 = @(t) t(ismember(t.Seed, 91001:91010), :);
Tsitl = table();
add = @(T, c, d, m, s) [T; table(string(c), string(d), m(1), m(2), s(1), s(2), 'VariableNames', {'Case','Design','MATLAB','MATLAB_SD','SITL','SITL_SD'})];
Tsitl = add(Tsitl, "S6", "PID-Budget", ms(C1.RMSE_mm(C1.Experiment=="E1")), ms(A.RMSE_mm(A.Group=="PBS")));
for d = ["AA-Budget","AA-HalfGain"], v = V(V.Exp=="E1" & V.Ctrl==d,:); Tsitl = add(Tsitl, "S6", d, [v.RMSE_ML v.SD_ML], [v.RMSE_SITL v.SD_SITL]); end
x = rd('matlab_ROx.csv'); Tsitl = add(Tsitl, "S6", "RO-HalfGain", ms(x.RMSE_mm), ms(A.RMSE_mm(A.Group=="ROS")));
x = rd('matlab_PIDS.csv'); Tsitl = add(Tsitl, "S6 strict", "PID-Budget-S", ms(x.RMSE_mm), ms(A.RMSE_mm(A.Group=="PSS")));
for d = ["AA-Budget-Strict","AA-HalfGain-Strict"], v = V(V.Exp=="E0" & V.Ctrl==d,:); Tsitl = add(Tsitl, "S6 strict", d, [v.RMSE_ML v.SD_ML], [v.RMSE_SITL v.SD_SITL]); end
b = C1r(C1r.Design=="PID-Budget" & C1r.Scen==5 & C1r.Config=="b_star",:); Tsitl = add(Tsitl, "S5", "PID-Budget", ms(b.RMSE_mm), ms(A.RMSE_mm(A.Group=="PB5")));
for d = ["AA-Budget","AA-HalfGain"], v = V(V.Exp=="S5" & V.Ctrl==d,:); Tsitl = add(Tsitl, "S5", d, [v.RMSE_ML v.SD_ML], [v.RMSE_SITL v.SD_SITL]); end
x = first10(rd('matlab_ROx_S5.csv')); Tsitl = add(Tsitl, "S5", "RO-HalfGain", ms(x.RMSE_mm), ms(A.RMSE_mm(A.Group=="RO5")));
x = first10(rd('matlab_A1x.csv')); Tsitl = add(Tsitl, "S6 ablation", "A1-Commanded", ms(x.RMSE_mm), ms(A.RMSE_mm(A.Group=="A1X")));
for d = ["A2-Unaligned","A3-NoIncrement"], v = V(V.Exp=="E3" & V.Ctrl==d,:); Tsitl = add(Tsitl, "S6 ablation", d, [v.RMSE_ML v.SD_ML], [v.RMSE_SITL v.SD_SITL]); end
Tsitl.AbsDiff = abs(Tsitl.SITL - Tsitl.MATLAB); writetable(Tsitl, fullfile(RA,'sitl_vs_matlab.csv'));
% ---- stack hover table (Ms per class with the 6-ms lag; hover predictions/measurements) ----
E = rd('designs_E_verify.csv'); N = rd('check_nominal_on_stack.csv'); I = rd('est_design_iteration0.csv'); H = rd('sitl_hover_E.csv'); HC = rd('increment_hover.csv');
Thover = table(); cls = ["xy","z","phth","psi"];
for d = ["PID-Budget","PID-E","AA-Budget-E","AA-HalfGain-E"]
    if d == "PID-Budget", src = N(N.Design=="PID-Budget",:); else, src = E(E.Design==d,:); end
    msE = arrayfun(@(c) src.Ms_E(src.Class==c), cls);
    if d == "PID-Budget", pc = src.sigma_du_est(src.Class=="xy"); pp = src.sigma_du_est(src.Class=="psi"); pj = NaN; mrp = NaN; mpsi = NaN; mz = NaN; mj = NaN;
    else, ii = I(I.Design==d,:); pc = ii.Pred_c_designtime(ii.Quantity=="sdu_rp"); pp = ii.Pred_c_designtime(ii.Quantity=="sdu_psi"); pj = ii.Pred_c_designtime(ii.Quantity=="J");
         h = H(H.Design==d,:); mrp = h.SdDU_rp; mpsi = h.SdDU_psi; mz = h.SdDU_z; mj = h.J; end
    hc = HC(HC.Design==d,:);
    Thover = [Thover; table(d, msE(1), msE(2), msE(3), msE(4), pc, hc.Pooled, hc.CI, pp, mpsi, mz, pj, mj, ...
        'VariableNames', {'Design','Ms_xy','Ms_z','Ms_phth','Ms_psi','PredCas','MeasCas','MeasCasCI','PredPsi','MeasPsi','MeasZ','PredJ','MeasJ'})]; %#ok<AGROW>
end
writetable(Thover, fullfile(RA,'stack_hover_table.csv'));
% ---- stack tracking table ----
S = rd('increment_sitl.csv'); Ttrack = S; writetable(Ttrack, fullfile(RA,'stack_tracking_table.csv'));
disp(Tsitl); disp(Thover);
end
