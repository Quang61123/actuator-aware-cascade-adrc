function clean_regeneration_test()
%CLEAN_REGENERATION_TEST Copies only the scripts and the result files needed for the SITL/stack summary tables
%   (sitl_vs_matlab, stack_hover_table, stack_tracking_table) and the budget-performance figure into an empty folder, regenerates them there (fresh MATLAB session, original folders not on
%   the path) and compares them with the files in the working tree: column names, row count and row keys,
%   text columns, NaN/Inf masks and the relative difference of the finite numeric entries.
%   Report: results/analysis/clean_regeneration_test.txt
src = fullfile(fileparts(mfilename('fullpath')), '..'); dst = fullfile(tempdir, ['clean_regen_' datestr(now,'yyyymmddHHMMSS')]);
mk = @(d) mkdir(fullfile(dst, d)); cp = @(f) copyfile(fullfile(src, f), fullfile(dst, f));
cellfun(mk, {'analysis','MATLAB_simulation','PX4_SITL\results','results\analysis\runs','results\estimator_noise\runs'});
cellfun(cp, {'analysis\regenerate_tables.m','analysis\sitl_increment_summary.m','analysis\make_figures_setpoint.m','MATLAB_simulation\sim_config.m', ...
    'PX4_SITL\results\sitl_vs_matlab_summary.csv','PX4_SITL\results\sitl_all_runs.csv','results\trials.csv'});
rf = {'setpoint_trials','setpoint_weight_runs','matlab_PIDS','matlab_ROx','matlab_ROx_S5','matlab_A1x','sitl_all_runs', ...
    'designs_E_verify','check_nominal_on_stack','est_design_iteration0','sitl_hover_E'};
for k = 1:numel(rf), cp(fullfile('results','analysis',[rf{k} '.csv'])); end
h = [dir(fullfile(src,'results','analysis','runs','E0H_hover_*.mat'))]; for k = 1:numel(h), copyfile(fullfile(h(k).folder,h(k).name), fullfile(dst,'results','analysis','runs',h(k).name)); end
h = dir(fullfile(src,'results','estimator_noise','runs','est_hover_c1_seed*.mat')); for k = 1:numel(h), copyfile(fullfile(h(k).folder,h(k).name), fullfile(dst,'results','estimator_noise','runs',h(k).name)); end
cmd = sprintf('matlab -batch "cd(''%s''); addpath(''%s'',''%s''); sitl_increment_summary; regenerate_tables(''%s''); make_figures_setpoint;" > "%s" 2>&1', ...
    fullfile(dst,'analysis'), fullfile(dst,'analysis'), fullfile(dst,'MATLAB_simulation'), dst, fullfile(dst,'clean.log'));
st = system(cmd);
out = {'sitl_vs_matlab.csv','stack_hover_table.csv','stack_tracking_table.csv','increment_sitl.csv','increment_hover.csv','budget_performance.csv'};
fid = fopen(fullfile(src,'results','analysis','clean_regeneration_test.txt'),'w'); [~, dname] = fileparts(dst); fprintf(fid,'clean folder: <tempdir>\\%s\nexit status: %d\n', dname, st);
nfail = 0;
for k = 1:numel(out)
    a = fullfile(dst,'results','analysis',out{k}); b = fullfile(src,'results','analysis',out{k});
    if ~exist(a,'file'), fprintf(fid,'%s: NOT REGENERATED\n', out{k}); nfail = nfail + 1; continue, end
    if ~exist(b,'file'), fprintf(fid,'%s: regenerated (no reference copy)\n', out{k}); nfail = nfail + 1; continue, end
    [ok, msg] = compare_tables(readtable(a,'TextType','string'), readtable(b,'TextType','string'));
    nfail = nfail + ~ok; fprintf(fid,'%s: %s\n', out{k}, msg);
end
fprintf(fid,'result: %s (%d of %d files differ)\n', string(ifelse(nfail == 0, "PASS", "FAIL")), nfail, numel(out));
fclose(fid); type(fullfile(src,'results','analysis','clean_regeneration_test.txt'));
end

function [ok, msg] = compare_tables(A, B)
% Compare a regenerated table A with the reference B.
if ~isequal(A.Properties.VariableNames, B.Properties.VariableNames), ok = false; msg = 'FAIL column names differ'; return, end
if height(A) ~= height(B), ok = false; msg = sprintf('FAIL row count %d vs %d', height(A), height(B)); return, end
isnum = varfun(@isnumeric, B, 'OutputFormat', 'uniform'); istxt = ~isnum;
% row keys and all text columns
txtok = true;
for j = find(istxt)
    x = string(A{:,j}); y = string(B{:,j}); x(ismissing(x)) = "<missing>"; y(ismissing(y)) = "<missing>";
    txtok = txtok && isequal(x, y);
end
% numeric columns: NaN/Inf masks, then finite values
na = A{:, isnum}; nb = B{:, isnum};
maskok = isequal(isnan(na), isnan(nb)) && isequal(isinf(na) & na > 0, isinf(nb) & nb > 0) && isequal(isinf(na) & na < 0, isinf(nb) & nb < 0);
fin = isfinite(na) & isfinite(nb);
d = 0; if any(fin(:)), d = max(abs(na(fin) - nb(fin)) ./ max(1, abs(nb(fin)))); end
ok = txtok && maskok && d <= 1e-12;
msg = sprintf('%s rows %d, text/key columns %s, NaN/Inf masks %s, max relative difference %.3g', ...
    char(ifelse(ok, "PASS", "FAIL")), height(B), char(ifelse(txtok, "equal", "DIFFER")), char(ifelse(maskok, "equal", "DIFFER")), d);
end

function y = ifelse(c, a, b)
if c, y = a; else, y = b; end
end
