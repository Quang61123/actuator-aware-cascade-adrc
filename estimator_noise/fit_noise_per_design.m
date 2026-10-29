function NF = fit_noise_per_design(dsig)
%FIT_NOISE_PER_DESIGN Estimator-error statistics and coloured-noise models.
%   dsig (1x12): signal delays [samples] from estimate_delays (median over the S6
%   runs). Data: noise-only hover SITL runs with EKF2 input (est_hover_*),
%   six designs, seeds 96001-96010, window 20-60 s, controller rate 1 kHz.
%   Error vector of each channel class (as read by the class model
%   (chan_model_est); the y-phi cascade is written as y'' = g*theta' with
%   theta' = -phi):
%     X : [e_x  e_vx  e_theta  e_thetadot]     Y : [e_y  e_vy  -e_phi  -e_phidot]
%     Z : [e_z  e_vz]                          P : [e_psi  e_psidot]
%   Statistics outputs (../results/estimator_noise/):
%     est_noise_by_ctrl.csv   variance of e and of its 1-ms increment per
%                            design (controller dependence)
%     est_noise_corr.csv      zero-lag correlation within each class
%     est_psd_by_ctrl.png, est_xcorr.png
%   Model fit: a multi-output (vector) autoregressive time-series model of order n,
%   A(q) n(k) = e(k) (least squares over all training segments), is fitted to the
%   outputs of each class, normalized by their standard deviation, on the
%   training seeds 96001-96005, and written in innovations state-space form
%       z(k+1) = Af z(k) + Bf v(k),  n(k) = Cf z(k) + Df v(k),  v ~ N(0,I),
%   with Bf = K*L, Df = L, L*L' = innovation covariance (state dimension p*n).
%   Order rule: n = 1..15 minimizing the Bayesian information criterion on
%   the training data among the stationary fits (spectral radius < 1). (A least-squares AR model
%   reproduces the lag-0 and lag-1 covariances, hence the variance and the
%   increment variance of the training data, for any order; the order
%   therefore governs the spectral shape, which is checked on the PSD and on
%   the held-out seeds.)
%   Design-dependence rule: models are fitted per design if, for any class
%   output, the design-to-design spread (max/min) of the variance or of the
%   increment variance exceeds 1.2; a model pooled over all designs is
%   fitted in any case (NF.<class>_pooled) for comparison.
%   Validation on the held-out seeds 96006-96010: variance and increment
%   variance errors (%), est_noise_validation.csv, est_psd_fit.png.
%   Saves NF (models, delays, rules, scales) to est_noise_models.mat.
here = fileparts(mfilename('fullpath'));
DW = fullfile(px4sitl_work_dir(),'est_data');
RES = fullfile(here,'..','results','estimator_noise'); if ~exist(RES,'dir'), mkdir(RES); end
Ts = 1e-3; win = [20 60];
cls = struct('name', {'X','Y','Z','P'}, 'col', {[1 4 8 11], [2 5 7 10], [3 6], [9 12]}, ...
    'sgn', {[1 1 1 1], [1 1 -1 -1], [1 1], [1 1]}, ...
    'lab', {{'e_x','e_vx','e_theta','e_thetadot'}, {'e_y','e_vy','-e_phi','-e_phidot'}, {'e_z','e_vz'}, {'e_psi','e_psidot'}});
ids = [1 2 4 8 9 10]; train = 96001:96005; valid = 96006:96010;
% ---------------- load all hover runs ----------------
f = dir(fullfile(DW,'est_hover_*.mat'));
Dat = struct('ctrl',{},'seed',{},'e',{},'seg',{});
for r = 1:numel(f)
    E = estimation_errors(fullfile(DW,f(r).name), dsig, win);
    Dat(end+1) = struct('ctrl', E.meta.ctrl, 'seed', E.meta.seed, 'e', E.e, 'seg', E.seg); %#ok<AGROW>
end
assert(~isempty(Dat), 'fit_noise_per_design: no hover runs');
% ---------------- statistics per design ----------------
Tb = table(); Tc = table();
for c = 1:numel(cls)
    for j = ids
        dd = Dat([Dat.ctrl] == j);
        if isempty(dd), continue; end
        [v, vd] = datastats(dd, cls(c)); p = numel(v);
        Tb = [Tb; table(repmat(string(cls(c).name),p,1), repmat(j,p,1), string(cls(c).lab(:)), v(:), vd(:), ...
            'VariableNames', {'Class','Ctrl','Output','Var','IncVar'})]; %#ok<AGROW>
    end
    Y = classdata(Dat, cls(c)); Yc = cat(1, Y{:});
    Rm = corrcoef(Yc);
    for a = 1:numel(cls(c).col), for b = a+1:numel(cls(c).col)
        Tc = [Tc; table(string(cls(c).name), string(cls(c).lab{a}), string(cls(c).lab{b}), Rm(a,b), ...
            'VariableNames', {'Class','A','B','Corr'})]; %#ok<AGROW>
    end, end
end
writetable(Tb, fullfile(RES,'est_noise_by_ctrl.csv')); writetable(Tc, fullfile(RES,'est_noise_corr.csv'));
% design dependence: spread max/min over designs of each output's variance and increment variance
dep = false(1,numel(cls)); spread = zeros(1,numel(cls));
for c = 1:numel(cls)
    for a = 1:numel(cls(c).lab)
        B = Tb(Tb.Class == cls(c).name & Tb.Output == cls(c).lab{a},:);
        spread(c) = max([spread(c), max(B.Var)/min(B.Var), max(B.IncVar)/min(B.IncVar)]);
    end
    dep(c) = spread(c) > 1.2;
end
% PSD per design (class X) and cross-correlation figure
fig = figure('Visible','off','Position',[50 50 1100 700]);
for a = 1:4
    subplot(2,2,a);
    for j = ids
        dd = Dat([Dat.ctrl] == j & ismember([Dat.seed], train)); if isempty(dd), continue; end
        Y = classdata(dd, cls(1)); [pp, fq] = segpsd(Y, a, Ts); loglog(fq, pp); hold on
    end
    grid on; title(cls(1).lab{a}, 'Interpreter','none'); xlabel('f [Hz]'); ylabel('PSD');
end
legend(arrayfun(@(j) sprintf('c%d',j), ids, 'uni', 0), 'Location','southwest');
exportgraphics(fig, fullfile(RES,'est_psd_by_ctrl.png'), 'Resolution', 130); close(fig);
fig = figure('Visible','off','Position',[50 50 900 350]);
Y = classdata(Dat(ismember([Dat.seed], train)), cls(1)); L = 200;
pairs = [1 2; 3 4; 1 3; 2 4]; for a = 1:4
    subplot(1,4,a); cc = segxcorr(Y, pairs(a,1), pairs(a,2), L); plot((-L:L)*Ts*1e3, cc); grid on
    title(sprintf('%s vs %s', cls(1).lab{pairs(a,1)}, cls(1).lab{pairs(a,2)}), 'Interpreter','none'); xlabel('lag [ms]');
end
exportgraphics(fig, fullfile(RES,'est_xcorr.png'), 'Resolution', 130); close(fig);
% ---------------- model fits ----------------
NF = struct('dsig', dsig, 'Ts', Ts, 'train', train, 'valid', valid, 'spread', spread, 'perdesign', dep, 'cls', cls);
Tv = table(); To = table(); lg = fullfile(RES,'est_fit_progress.txt'); fid = fopen(lg,'w'); fclose(fid);
for c = 1:numel(cls)
    groups = {ids}; if dep(c), groups = [{ids}, num2cell(ids)]; end    % pooled model first
    for g = 1:numel(groups)
        dtr = Dat(ismember([Dat.ctrl], groups{g}) & ismember([Dat.seed], train));
        dva = Dat(ismember([Dat.ctrl], groups{g}) & ismember([Dat.seed], valid));
        Ytr = classdata(dtr, cls(c)); Yc = cat(1, Ytr{:}); p = size(Yc,2);
        sc = std(Yc, 1, 1); Ytrn = cellfun(@(y) (y - mean(y,1))./sc, Ytr, 'uni', 0);
        [vtr, vdtr] = datastats(dtr, cls(c));
        bic = inf(1,15); mods = cell(1,15);
        for n = 1:15
            [mods{n}, bic(n)] = varfit(Ytrn, n);
            To = [To; table(string(cls(c).name), string(mat2str(groups{g})), n, bic(n), 'VariableNames', {'Class','Designs','Order','BIC'})]; %#ok<AGROW>
        end
        % lowest BIC among the stationary fits (spectral radius of Af < 1)
        [~, ord] = sort(bic); chosen = NaN;
        for n = ord
            cand = var2filter(mods{n}, sc);
            if max(abs(eig(cand.Af))) < 1, chosen = n; best = cand; break; end
        end
        assert(~isnan(chosen), 'fit_noise_per_design: no stationary VAR fit for class %s', cls(c).name);
        best.order = chosen; best.ctrls = groups{g}; best.scale = sc;
        [vm, vdm] = filterstats(best); [vva, vdva] = datastats(dva, cls(c));
        Tv = [Tv; table(repmat(string(cls(c).name),p,1), repmat(string(mat2str(groups{g})),p,1), repmat(chosen,p,1), string(cls(c).lab(:)), ...
            100*(vm(:)./vtr(:) - 1), 100*(vdm(:)./vdtr(:) - 1), 100*(vm(:)./vva(:) - 1), 100*(vdm(:)./vdva(:) - 1), ...
            'VariableNames', {'Class','Designs','Order','Output','TrainVarErrPct','TrainIncVarErrPct','ValidVarErrPct','ValidIncVarErrPct'})]; %#ok<AGROW>
        if g == 1, NF.([cls(c).name '_pooled']) = best; end
        if ~dep(c) || g > 1
            if ~isfield(NF, cls(c).name), NF.(cls(c).name) = best; else, NF.(cls(c).name)(end+1) = best; end
        end
        fid = fopen(lg,'a'); fprintf(fid,'%s class %s designs %s order %d\n', char(datetime('now','Format','HH:mm:ss')), cls(c).name, mat2str(groups{g}), chosen); fclose(fid);
    end
end
writetable(Tv, fullfile(RES,'est_noise_validation.csv')); writetable(To, fullfile(RES,'est_noise_order_search.csv'));
% PSD fit figure (validation seeds)
fig = figure('Visible','off','Position',[50 50 1200 800]); pn = 0;
for c = 1:numel(cls)
    m = NF.(cls(c).name)(1); dva = Dat(ismember([Dat.ctrl], m.ctrls) & ismember([Dat.seed], valid));
    Y = classdata(dva, cls(c)); sysm = ss(m.Af, m.Bf, m.Cf, m.Df, Ts);
    for a = 1:numel(cls(c).col)
        pn = pn + 1; subplot(3,4,pn); [pp, fq] = segpsd(Y, a, Ts); loglog(fq, pp, 'Color', [.6 .6 .6]); hold on
        H = squeeze(freqresp(sysm(a,:), 2*pi*fq)); pm = 2*Ts*sum(abs(reshape(H, [], numel(fq))).^2, 1);
        loglog(fq, pm, 'k', 'LineWidth', 1.2); grid on; title(sprintf('%s: %s', cls(c).name, cls(c).lab{a}), 'Interpreter','none');
    end
end
exportgraphics(fig, fullfile(RES,'est_psd_fit.png'), 'Resolution', 130); close(fig);
save(fullfile(RES,'est_noise_models.mat'), 'NF');
end

% ------------------------------------------------------------------------
function Y = classdata(dd, cl)
% cell array of contiguous segments (N x p) of the class error vector
Y = {};
for r = 1:numel(dd)
    e = dd(r).e(:, cl.col).*cl.sgn;
    for s = unique(dd(r).seg)'
        y = e(dd(r).seg == s, :);
        if size(y,1) > 1000, Y{end+1} = y; end %#ok<AGROW>
    end
end
end

function [v, vd] = datastats(dd, cl)
% pooled variance of the class outputs and of their 1-ms increments (per-run mean removed)
Y = classdata(dd, cl); n = 0; s = 0; nd = 0; sd = 0;
for i = 1:numel(Y)
    y = Y{i} - mean(Y{i},1); d = diff(Y{i});
    s = s + sum(y.^2,1); n = n + size(y,1); sd = sd + sum(d.^2,1); nd = nd + size(d,1);
end
v = s/n; vd = sd/nd;
end

function [V, bic] = varfit(Y, n)
% least-squares vector AR(n) fit, n(k) = A1 n(k-1) + ... + An n(k-n) + e(k),
% over all segments (each segment contributes its own lagged regressions)
p = size(Y{1},2); G = zeros(n*p); H = zeros(n*p, p); N = 0;
for i = 1:numel(Y)
    y = Y{i}; m = size(y,1);
    Z = zeros(m-n, n*p);
    for l = 1:n, Z(:, (l-1)*p+(1:p)) = y(n+1-l:m-l, :); end
    G = G + Z'*Z; H = H + Z'*y(n+1:m,:); N = N + m - n;
end
B = G\H;                                   % (n*p) x p, rows = lag blocks
Sig = 0;
for i = 1:numel(Y)
    y = Y{i}; m = size(y,1); Z = zeros(m-n, n*p);
    for l = 1:n, Z(:, (l-1)*p+(1:p)) = y(n+1-l:m-l, :); end
    E = y(n+1:m,:) - Z*B; Sig = Sig + E'*E;
end
Sig = Sig/N;
V = struct('A', B', 'Sig', Sig, 'n', n, 'p', p);   % A = [A1 ... An] (p x n*p)
bic = N*log(det(Sig)) + n*p^2*log(N);
end

function m = var2filter(V, sc)
% VAR model of normalized outputs -> innovations state space -> shaping filter
p = V.p; n = V.n; Lm = chol(V.Sig, 'lower'); S = diag(sc);
Af = [V.A; eye(p*(n-1)), zeros(p*(n-1), p)];   % state z(k) = [n(k-1); ...; n(k-n)]
Bf = [Lm; zeros(p*(n-1), p)];
m = struct('Af', Af, 'Bf', Bf, 'Cf', S*V.A, 'Df', S*Lm, 'order', NaN, 'ctrls', [], 'scale', sc);
end

function [v, vd] = filterstats(m)
% stationary variance of n and of its 1-ms increment
P = dlyap(m.Af, m.Bf*m.Bf');
V = m.Cf*P*m.Cf' + m.Df*m.Df';
H = m.Cf*(m.Af - eye(size(m.Af,1))); Q = m.Cf*m.Bf - m.Df;
Vd = H*P*H' + Q*Q' + m.Df*m.Df';
v = diag(V)'; vd = diag(Vd)';
end

function [pp, fq] = segpsd(Y, a, Ts)
nf = 4096; acc = 0; nseg = 0;
for i = 1:numel(Y)
    y = Y{i}(:,a); if numel(y) < nf, continue; end
    [p, fq] = pwelch(y - mean(y), nf, [], [], 1/Ts); acc = acc + p*numel(y); nseg = nseg + numel(y);
end
pp = acc/nseg; fq = fq(:); pp = pp(:);
pp = pp(2:end); fq = fq(2:end);
end

function cc = segxcorr(Y, a, b, L)
acc = 0; nn = 0;
for i = 1:numel(Y)
    x = Y{i}(:,a) - mean(Y{i}(:,a)); y = Y{i}(:,b) - mean(Y{i}(:,b));
    c = xcorr(x, y, L, 'coeff'); acc = acc + c*numel(x); nn = nn + numel(x);
end
cc = acc/nn;
end
