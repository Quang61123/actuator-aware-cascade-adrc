function M = results_manifest()
%RESULTS_MANIFEST Result item -> script -> result file -> design/gain source -> SHA-256 of the
%   result file. Writes results/analysis/results_manifest.csv.
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
L = {
 'Gains: strict PID, -E and b* designs', 'design_pid_strict.m, design_estimator_stack.m, finalize_estimator_designs.m, setpoint_weight_select.m', 'results/analysis/pidS.mat', 'design_gains.mat C(1:7), Strict; pidS.mat; designs_E.mat; setpoint_designs.mat';
 'Gains: estimator-stack (-E) designs', 'finalize_estimator_designs.m', 'results/analysis/designs_E.mat', 'designs_E.mat DE(1:3)';
 'PID set-point weights b*', 'setpoint_weight_select.m (rule setpoint_rule.txt)', 'results/analysis/setpoint_designs.mat', 'setpoint_designs.mat C(k).star';
 'Response requirement check', 'settling_check.m', 'results/analysis/settling_check.csv', 'all designs';
 'MATLAB trials, PID-Budget (b*)', 'setpoint_weight_trials.m', 'results/analysis/setpoint_trials.csv', 'PID-Budget (b*)';
 'MATLAB trials, primary designs', 'MATLAB_simulation/mc_trials.m', 'results/trials.csv', 'design_gains.mat';
 'MATLAB trials, RO-HalfGain', 'matlab_design_runs.m (ROx), matlab_extra_runs.m', 'results/analysis/matlab_ROx.csv', 'designs_main.mat M.ROx';
 'MATLAB trials, PID-Budget-S', 'matlab_design_runs.m (PIDS)', 'results/analysis/matlab_PIDS.csv', 'pidS.mat';
 'Ablation A1-Commanded', 'matlab_design_runs.m (A1x), matlab_paired.m', 'results/analysis/matlab_A1x.csv', 'designs_main.mat M.A1x';
 'Ablation A4-CommandedTilt', 'tilt_input_ablation.m, tilt_input_scan.m', 'results/analysis/tilt_scan_summary.csv', 'tilt_designs.mat; tilt_scan_verify.csv';
 'Commanded-input ratio sensitivity', 'commanded_low_ratio.m', 'results/analysis/lowratio_summary.csv', 'commanded_scan_ratio068.mat';
 'Budget-performance curve', 'make_figures_setpoint.m', 'results/analysis/budget_performance.csv', 'trials E4 (AA); setpoint_trials E4 (PID b*)';
 'Noise-activity and mismatch-count figures', 'make_figures_setpoint.m', 'results/analysis/fig_mismatch_counts.pdf', 'trials E5/E7 (AA); setpoint_trials (PID b*)';
 'SITL vs MATLAB summary', 'regenerate_tables.m', 'results/analysis/sitl_vs_matlab.csv', 'PX4 ids 2,4,9,10 (first campaign), 15, 20, 21, 22';
 'SITL lockstep misses', 'lockstep_miss_sensitivity.m, report_sitl_commanded.m, sitl_increment_summary.m', 'results/analysis/miss_by_group.csv', 'run plans run_plan*.csv';
 'Noise-model validation (3 levels)', 'validate_noise_models.m, noise_delay_consistency.m', 'results/analysis/validation_3level.csv', 'noise_long/, noise_est_data_rateconv/';
 'Stack hover table', 'regenerate_tables.m (check_nominal_on_stack.m, est_design_iteration.m, sitl_increment_summary.m)', 'results/analysis/stack_hover_table.csv', 'designs_E.mat';
 'Stack tracking table', 'sitl_increment_summary.m', 'results/analysis/increment_sitl.csv', 'PX4 ids 17-19, 22';
 'Frozen ceilings', 'design_estimator_stack.m', 'results/analysis/frozen_ceilings.csv', 'PID-Budget, pooled estimator-campaign model';
 'PID equation-code consistency', 'analysis/pid_equation_check.m', 'results/analysis/pid_equation_consistency.csv', 'setpoint_designs.mat C(k).star';
};
M = cell2table(L, 'VariableNames', {'Item','Script','ResultFile','DesignSource'}); M.SHA256 = strings(height(M),1); M.Bytes = nan(height(M),1);
for k = 1:height(M)
    f = fullfile(root, strrep(M.ResultFile{k}, '/', filesep)); if ~exist(f,'file'), M.SHA256(k) = "MISSING"; continue, end
    fid = fopen(f,'r'); b = fread(fid, inf, '*uint8'); fclose(fid); md = java.security.MessageDigest.getInstance('SHA-256'); md.update(b);
    h = typecast(md.digest(), 'uint8'); M.SHA256(k) = string(lower(reshape(dec2hex(h)', 1, []))); M.Bytes(k) = numel(b);
end
writetable(M, fullfile(root,'results','analysis','results_manifest.csv')); disp(M(:, {'Item','ResultFile','Bytes'}));
end
