function DE = design_estimator_initial()
%DESIGN_ESTIMATOR_INITIAL Equal-activity designs for the estimator-fed loop.
%   Ceiling: capE = estimator-aware predicted cascade inner increment of
%   PID-Budget (noise models pooled over the designs, est_stats), i.e. the
%   equal-activity rule of the primary designs applied to the estimator-fed loop.
%   AA-Budget and AA-HalfGain are redesigned under capE (names -E); PID-Budget
%   is unchanged. The outer x-y loop alone cannot meet capE (the inner
%   increment is dominated by the filtered-rate noise passing through the
%   inner loop), so the inner roll/pitch pair (omega_c, omega_o) is searched
%   as well, on the grid omega_c = logspace(0, log10(omega_c,nominal), 10),
%   omega_o/omega_c = logspace(log10(1.5), log10(60), 12), omega_o <= S.wo_max:
%     (i)  the inner roll/pitch class must meet its local ceilings of the
%          primary design (white model), Ms <= S.Msmax and OS <= S.osmax;
%     (ii) a prefilter discards inner pairs whose cascade increment already
%          exceeds capE at the lowest outer bandwidth of the search grid;
%     (iii) cascade_outer_design_E then selects the outer loop under capE
%          (all other outer constraints as in the primary design).
%   The design with the smallest estimator-aware tracking standard deviation
%   is kept. Saves ../results/estimator_noise/est_designs_E.mat (DE, capE,
%   ledger, info, search) and est_designs_E.csv, est_design_E_search.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RES = fullfile(RP,'estimator_noise');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub;
L = load(fullfile(RES,'est_noise_models.mat')); NF = L.NF;
zp = est_stats(D.C(1), S, NF); capE = zp(4);
lg = fullfile(RES,'est_design_E_progress.txt'); fid = fopen(lg,'w'); fprintf(fid,'capE = %.6g\n', capE); fclose(fid);
DE = D.C([2 4]); ledger = table(); info = table(); search = table();
for j = 1:2
    base = DE(j); bestCost = inf; bestCo = [];
    wcs = logspace(0, log10(base.G(19)), 10); rts = logspace(log10(1.5), log10(60), 12);
    for wc = wcs
        for rt = rts
            wo = rt*wc; if wo > S.wo_max, continue, end
            c = base; gi = gains_of(c.kind, rt, wc)';
            c.G(19:24) = gi; c.G(25:30) = gi;
            g = c.G(19:24)'; oa = [c.tauhat c.ndhat];
            [~,~,vu,vdu] = chan_var(c.kind, g, S.Ts, S.sy(3), S.sv(3), 0, [], 0, 0, S.act, oa);
            msi = chan_ms(c.kind, g, S.Ts, [], S.act, [], oa, 600);
            okIn = isfinite(vu) && sqrt(vu) <= S.ub(3)*(1+1e-9) && sqrt(vdu) <= dub(3)*(1+1e-9) && msi <= S.Msmax;
            if okIn, osi = chan_step(c.kind, g, S.Ts, [], S.act, 4*S.Tstep(3), oa); okIn = osi <= S.osmax; end
            if ~okIn, continue, end
            c0 = c; c0.G(1:6) = gains_of(c0.kind, 10, 0.05)'; c0.G(7:12) = c0.G(1:6); z0 = est_stats(c0, S, NF);
            if z0(4) > capE*(1-0.01), continue, end
            c.name = [base.name '-E'];
            [co, Lg] = cascade_outer_design_E(c, S, dub, capE, NF);
            if isempty(co), search = [search; table(string(c.name), wc, wo, NaN, NaN, NaN, 'VariableNames', {'Design','Inner_wc','Inner_wo','Outer_wc','SigmaE_est','SigmaDU_est'})]; continue, end %#ok<AGROW>
            ze = est_stats(co, S, NF);
            search = [search; table(string(c.name), wc, wo, co.G(1), ze(1), ze(4), 'VariableNames', {'Design','Inner_wc','Inner_wo','Outer_wc','SigmaE_est','SigmaDU_est'})]; %#ok<AGROW>
            logline(lg, sprintf('%s inner wc %.3f wo %.3f -> outer wc %.4f, sigmaE_est %.4f, sigmaDU_est %.5f', c.name, wc, wo, co.G(1), ze(1), ze(4)));
            if ze(1) < bestCost, bestCost = ze(1); bestCo = co; bestL = Lg; end
        end
    end
    if isempty(bestCo), logline(lg, sprintf('%s: NO FEASIBLE -E DESIGN', base.name)); DE(j).name = [base.name '-E-INFEASIBLE']; continue, end
    co = bestCo; ze = est_stats(co, S, NF); zw = cascade_stats(co, S);
    msc = chan_ms(co.kind, co.G(1:6)', S.Ts, co.G(19:24)', S.act, [], [co.tauhat co.ndhat], 600);
    msi = chan_ms(co.kind, co.G(19:24)', S.Ts, [], S.act, [], [co.tauhat co.ndhat], 600);
    logline(lg, sprintf('SELECTED %s: inner wc %.4g wo %.4g, outer wc %.4g wo %.4g, sigmaE_est %.4g, sigmaDU_est %.5g', ...
        co.name, co.G(19), co.G(20), co.G(1), co.G(2), ze(1), ze(4)));
    DE(j) = co; ledger = [ledger; bestL]; %#ok<AGROW>
    info = [info; table(string(co.name), co.G(1), co.G(2), co.G(19), co.G(20), ze(1), ze(4), capE, zw(1), zw(4), msc, msi, ...
        'VariableNames', {'Design','Outer_wc','Outer_wo','Inner_wc','Inner_wo','SigmaE_est','SigmaDU_est','CapE', ...
        'SigmaE_white','SigmaDU_white','Ms_cascade','Ms_inner'})]; %#ok<AGROW>
end
save(fullfile(RES,'est_designs_E.mat'), 'DE', 'capE', 'ledger', 'info', 'search');
writetable(info, fullfile(RES,'est_designs_E.csv')); writetable(search, fullfile(RES,'est_design_E_search.csv'));
end

function logline(f, s)
fid = fopen(f,'a'); fprintf(fid,'%s %s\n', char(datetime('now','Format','HH:mm:ss')), s); fclose(fid);
end
