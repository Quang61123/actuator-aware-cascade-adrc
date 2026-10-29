function make_figures_setpoint()
%MAKE_FIGURES_SETPOINT Figures with the b* set-point weights: budget-performance curve (AA points from trials E4, PID points from
%   setpoint_trials E4 with the b* weights), noise activity (AA from trials E5-1..4, PID-Budget (b*)
%   from setpoint_trials) and mismatch counts ( AA from trials E7, PID-Budget (b*) from setpoint_trials;
%   each cell annotated F/E = failures/bounded threshold exceedances). Style as make_figures.
%   Output: results/analysis/fig_*.pdf/png and budget_performance.csv.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); P = sim_config(); set(groot,'defaultFigureVisible','off');
O = rd(fullfile(RP,'trials.csv')); N = rd(fullfile(RA,'setpoint_trials.csv'));
% ---- budget-performance curve ----
F = table(); fam = {'PID','P2','HG-AA'}; lab = {'PID','AA-Budget','AA-HalfGain'}; gm = [.5 .75 1 1.25 1.5];
for k = 1:3, for g = gm
  if k == 1, a = N(N.Experiment=="E4" & abs(N.Gamma-g)<1e-9 & N.Config=="b_star",:); else, a = O(O.Experiment=="E4" & O.Ctrl==string(sprintf('%s-g%03d',fam{k},round(100*g))),:); end
  v = a.RMSE_mm(~a.Failed); n = numel(v); F = [F; table(string(lab{k}), g, n, mean(v), tinv(.975,n-1)*std(v)/sqrt(n), mean(a.Increment(~a.Failed)), 'VariableNames', {'Family','Gamma','N','RMSE','CI','Inc'})]; %#ok<AGROW>
end, end
writetable(F, fullfile(RA,'budget_performance.csv')); disp(F)
f = figure('Position',[100 100 800 340]); tiledlayout(1,2,'TileSpacing','compact'); col = lines(3);
for panel = 1:2, nexttile; hold on; for k = 1:3, z = F(F.Family==string(lab{k}),:); x = z.Gamma; if panel == 2, x = z.Inc; end
  errorbar(x, z.RMSE, z.CI, '-o', 'Color', col(k,:), 'LineWidth', 1, 'MarkerSize', 4); end
  ylabel('Position RMSE [mm]'); grid on; if panel == 1, xlabel('Local increment-ceiling multiplier \gamma (cascade ceiling fixed)'); else, xlabel('RMS rotor-command increment [rad/s]'); end; if panel == 1, lloc = 'northeast'; else, lloc = 'east'; end; legend(lab, 'Location', lloc, 'FontSize', 10); end
ex(f, RA, 'fig_budget_performance');
% ---- noise activity ----
f = figure('Position',[100 100 820 330]); tiledlayout(1,2,'TileSpacing','compact'); sc = [1 .5 2 4];
for panel = 1:2, nexttile; hold on; for c = {'PID-Budget','AA-Budget','AA-HalfGain'}
  yy = zeros(1,4); for z = 1:4, e = string(sprintf('E5-%d',z)); if strcmp(c{1},'PID-Budget'), a = N(N.Experiment==e,:); else, a = O(O.Experiment==e & O.Ctrl==string(c{1}),:); end
    if panel == 1, yy(z) = mean(a.RMSE_mm(~a.Failed)); else, yy(z) = mean(a.Increment(~a.Failed)); end, end
  [xs, ii] = sort(sc); plot(xs, yy(ii), '-o', 'DisplayName', c{1}); end
  grid on; xlabel('White sensor-noise standard-deviation multiplier'); if panel == 1, ylabel('Mean position RMSE [mm]'); else, ylabel('Rotor-command increment [rad/s]'); end; legend('Location','northwest','FontSize',9); end
ex(f, RA, 'fig_noise_activity');
% ---- mismatch counts ----
f = figure('Position',[100 100 980 300]); tiledlayout(1,3,'TileSpacing','compact'); names = {'PID-Budget','AA-Budget','AA-HalfGain'};
for c = 1:3, nexttile; X = nan(8,3); FF = X; EE = X;
  for it = 1:8, for id = 1:3
    if c == 1, a = N(N.Experiment=="E7" & abs(N.Tau-P.taus(it))<1e-9 & N.Delay==P.delays(id),:); else, a = O(O.Experiment=="E7" & O.Ctrl==string(names{c}) & abs(O.Tau-P.taus(it))<1e-9 & O.Delay==P.delays(id),:); end
    FF(it,id) = sum(a.Failed); EE(it,id) = sum(a.FluctuationThreshold); X(it,id) = FF(it,id) + EE(it,id); end, end
  imagesc(P.delays, P.taus*1000, X); set(gca,'YDir','normal'); clim([0 10]); colorbar; 
  for it = 1:8, for id = 1:3, if X(it,id) > 0, tc = 'k'; if X(it,id) <= 5, tc = 'w'; end; text(P.delays(id), P.taus(it)*1000, sprintf('%d/%d', FF(it,id), EE(it,id)), 'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',tc); end, end, end
  xlabel('True measurement delay [ms]'); ylabel('Rotor time constant [ms]'); title(names{c}); end
ex(f, RA, 'fig_mismatch_counts');
end
function T = rd(fn), T = readtable(fn,'TextType','string'); for v = {'Failed','FluctuationThreshold'}, x = T.(v{1}); if isstring(x) || iscell(x), T.(v{1}) = strcmpi(string(x),"true") | string(x)=="1"; end, end, if ~ismember('Config', T.Properties.VariableNames), T.Config = repmat("", height(T), 1); end, end
function ex(f, d, nm), exportgraphics(f, fullfile(d,[nm '.pdf']), 'ContentType','vector'); exportgraphics(f, fullfile(d,[nm '.png']), 'Resolution',180); close(f); end