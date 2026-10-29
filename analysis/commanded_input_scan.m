function [T, F] = commanded_input_scan(famSel, rg, tag)
%   Optional: famSel = row indices of the families,
%   rg = ratio grid, tag = suffix of the output files (sensitivity run).
if nargin < 1 || isempty(famSel), famSel = [1 2]; end
if nargin < 2 || isempty(rg), rg = logspace(log10(1.5), log10(60), 20); end
if nargin < 3, tag = ''; end
%COMMANDED_INPUT_SCAN Re-search of the commanded-input designs under the
%   response requirement (RO-HalfGain and A1-Commanded are searched again
%   under the settling constraint).
%   For each design family (A1-Commanded: kind 4; RO-HalfGain: kind 8; same
%   internal model, budgets and search ranges as the original designs) and
%   each inner class (z, phi/theta, psi), a dense grid over the observer/
%   controller ratio r = wo/wc (20 log points, 1.5-60, as the original
%   search) and wc (60 log points, 0.05 to wo_max/r) is evaluated against
%   ALL constraints of the covariance design problem plus the response requirement:
%     sigma_u <= ub(c), sigma_du <= dub(c) (local), Ms <= 1.98 (600-point
%     grid plus refinement in the search; candidates re-checked on 4000),
%     OS <= 10 % and 2% settling time <= req(c), both on a 400-s horizon.
%   The horizontal class is searched (grid as above, cascade cap 1.382 and
%   the same constraints) only when a feasible phi/theta design exists.
%   Outputs ../results/analysis/commanded_scan_grid.csv (all grid points)
%   and commanded_scan<tag>.csv/.mat (per family/class: feasible points, best
%   design, and the largest Ms-feasible wc with its settling time).
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub; cap = D.cap;
req = [15 5 2 5]; fam = {'A1-Commanded', 4; 'RO-HalfGain', 8}; fam = fam(famSel,:);
nw = 60;
F = table(); T = table();
for f = 1:size(fam,1)
    kind = fam{f,2}; oact = S.act; best = cell(1,4);
    for c = [2 3 4 1]
        if c == 1
            if isempty(best{3}), T = [T; summ(fam{f,1}, c, F(false(height(F),1),:), [])]; continue, end %#ok<AGROW>
            gin = best{3};
        else, gin = []; end
        rows = cell(numel(rg),1);
        parfor q = 1:numel(rg)
            r = rg(q); wg = logspace(log10(0.05), log10(S.wo_max/r), nw); R = zeros(nw, 11);
            for i = 1:nw
                g = gains_of(kind, r, wg(i));
                z = channel_stats(kind, g, c, S, gin, oact, S.act, false, 600);
                M = chan_model(kind, g, S.Ts, gin, S.act, [], oact); [os, ts] = step_response_metrics(M, S.Ts);
                cdu = NaN; okc = true;
                if c == 1
                    co = D.C(5); co.kind = kind; co.G(1:6) = g'; co.G(7:12) = g'; co.G(19:24) = gin'; co.G(25:30) = gin';
                    zz = cascade_stats(co, S); cdu = zz(4); okc = cdu <= cap*(1-0.01);
                end
                ok = all(isfinite(z(1:4))) && z(2) <= S.ub(c)*(1+1e-9) && z(3) <= dub(c)*(1+1e-9) && z(4) <= S.Msmax && z(6) < 1 ...
                    && os <= S.osmax && ts <= req(c) && okc;
                R(i,:) = [r, wg(i), r*wg(i), z(1), z(2), z(3), z(4), os, ts, cdu, ok];
            end
            rows{q} = R;
        end
        R = vertcat(rows{:});
        G = array2table(R, 'VariableNames', {'ratio','wc','wo','sigma_e','sigma_u','sigma_du','Ms','OS','ts','cascade_du','Feasible'});
        G.Design = repmat(string(fam{f,1}), height(G), 1); G.Class = repmat(c, height(G), 1);
        F = [F; G]; %#ok<AGROW>
        fe = G(G.Feasible == 1, :);
        if ~isempty(fe)
            [~, j] = min(fe.sigma_e); g = gains_of(kind, fe.ratio(j), fe.wc(j));
            ms4 = chan_ms(kind, g, S.Ts, gin, S.act, [], oact, 4000);
            if ms4 > S.Msmax + 1e-9  % re-check on the fine grid; drop failing points
                fe = fe(arrayfun(@(k) chan_ms(kind, gains_of(kind, fe.ratio(k), fe.wc(k)), S.Ts, gin, S.act, [], oact, 4000) <= S.Msmax, 1:height(fe)), :);
                if ~isempty(fe), [~, j] = min(fe.sigma_e); g = gains_of(kind, fe.ratio(j), fe.wc(j)); end
            end
            if ~isempty(fe), best{c} = g; end
        end
        T = [T; summ(fam{f,1}, c, G, best{c})]; %#ok<AGROW>
        fprintf('SCAN %s class %d: %d feasible of %d\n', fam{f,1}, c, sum(G.Feasible), height(G));
    end
end
writetable(F, fullfile(RA,['commanded_scan_grid' tag '.csv']));
writetable(T, fullfile(RA,['commanded_scan' tag '.csv']));
save(fullfile(RA,['commanded_scan' tag '.mat']), 'T', 'F', 'rg');
end

function row = summ(name, c, G, g)
% per family/class summary: feasibility and the binding constraint
nf = sum(G.Feasible == 1);
okL = G.sigma_u <= Inf; % placeholder for clarity
msok = G.Ms <= 1.98 & isfinite(G.Ms);
if any(msok)
    [~, j] = max(G.wc .* msok); tsMs = G.ts(j); wcMs = G.wc(j); rMs = G.ratio(j);
    tsBest = min(G.ts(msok));
else, tsMs = NaN; wcMs = NaN; rMs = NaN; tsBest = NaN; end
if isempty(g), gwc = NaN; gwo = NaN; else, gwc = g(1); gwo = g(2); end
row = table(string(name), c, height(G), nf, gwc, gwo, rMs, wcMs, tsMs, tsBest, ...
    'VariableNames', {'Design','Class','GridPoints','Feasible','Best_wc','Best_wo','Ratio_at_maxMsfeasible_wc','Max_Msfeasible_wc','ts_at_that_wc','Min_ts_with_Ms_ok'});
end
