function NF = fit_noise_models(files, ds, win, train, valid, outdir)
%FIT_NOISE_MODELS Estimator-error models for the stack designs, same method as
%   fit_noise_per_design (vector AR by least squares on normalized class outputs,
%   order 1-15 by BIC among stationary fits, innovations shaping filter),
%   generalized to any set of SITL hover logs:
%     files  : cell array of est_data/redesign_data .mat files (R, Lx, c, seed)
%     ds     : 1x12 delay convention [samples] (est_stack_tools('ds'))
%     win    : window [s]; train/valid: seed lists (disjoint)
%     outdir : output folder (NF in noise_models.mat, CSV tables)
%   Models: one per design (level a) and one pooled over all designs in the
%   files (level b), for every class X, Y, Z, P. Validation: variance and
%   increment-variance errors on the training and held-out seeds
%   (noise_validation.csv).
if ~exist(outdir,'dir'), mkdir(outdir); end
Ts = 1e-3;
cls = struct('name', {'X','Y','Z','P'}, 'col', {[1 4 8 11], [2 5 7 10], [3 6], [9 12]}, ...
    'sgn', {[1 1 1 1], [1 1 -1 -1], [1 1], [1 1]}, ...
    'lab', {{'e_x','e_vx','e_theta','e_thetadot'}, {'e_y','e_vy','-e_phi','-e_phidot'}, {'e_z','e_vz'}, {'e_psi','e_psidot'}});
Dat = struct('ctrl',{},'seed',{},'e',{},'seg',{});
for r = 1:numel(files)
    E = errs(files{r}, ds, win);
    Dat(end+1) = struct('ctrl', E.ctrl, 'seed', E.seed, 'e', E.e, 'seg', E.seg); %#ok<AGROW>
end
ids = unique([Dat.ctrl]);
NF = struct('dsig', ds, 'Ts', Ts, 'train', train, 'valid', valid, 'win', win, 'ids', ids, 'cls', cls);
Tv = table(); Tb = table();
for c = 1:numel(cls)
    for j = ids
        [v, vd] = datastats(Dat([Dat.ctrl] == j), cls(c)); p = numel(v);
        Tb = [Tb; table(repmat(string(cls(c).name),p,1), repmat(j,p,1), string(cls(c).lab(:)), v(:), vd(:), ...
            'VariableNames', {'Class','Ctrl','Output','Var','IncVar'})]; %#ok<AGROW>
    end
    groups = [{ids}, num2cell(ids)];
    for g = 1:numel(groups)
        dtr = Dat(ismember([Dat.ctrl], groups{g}) & ismember([Dat.seed], train));
        dva = Dat(ismember([Dat.ctrl], groups{g}) & ismember([Dat.seed], valid));
        Ytr = classdata(dtr, cls(c)); Yc = cat(1, Ytr{:}); p = size(Yc,2);
        sc = std(Yc, 1, 1); Ytrn = cellfun(@(y) (y - mean(y,1))./sc, Ytr, 'uni', 0);
        [vtr, vdtr] = datastats(dtr, cls(c));
        bic = inf(1,15); mods = cell(1,15);
        for n = 1:15, [mods{n}, bic(n)] = varfit(Ytrn, n); end
        [~, ord] = sort(bic); chosen = NaN;
        for n = ord
            cand = var2filter(mods{n}, sc);
            if max(abs(eig(cand.Af))) < 1, chosen = n; best = cand; break; end
        end
        assert(~isnan(chosen), 'fit_noise_models: no stationary VAR fit for class %s', cls(c).name);
        best.order = chosen; best.ctrls = groups{g}; best.scale = sc;
        [vm, vdm] = filterstats(best);
        if isempty(dva), vva = nan(1,p); vdva = vva; else, [vva, vdva] = datastats(dva, cls(c)); end
        Tv = [Tv; table(repmat(string(cls(c).name),p,1), repmat(string(mat2str(groups{g})),p,1), repmat(chosen,p,1), string(cls(c).lab(:)), ...
            100*(vm(:)./vtr(:) - 1), 100*(vdm(:)./vdtr(:) - 1), 100*(vm(:)./vva(:) - 1), 100*(vdm(:)./vdva(:) - 1), ...
            'VariableNames', {'Class','Designs','Order','Output','TrainVarErrPct','TrainIncVarErrPct','ValidVarErrPct','ValidIncVarErrPct'})]; %#ok<AGROW>
        if g == 1, NF.([cls(c).name '_pooled']) = best;
        elseif ~isfield(NF, cls(c).name), NF.(cls(c).name) = best;
        else, NF.(cls(c).name)(end+1) = best; end
    end
end
writetable(Tv, fullfile(outdir,'noise_validation.csv')); writetable(Tb, fullfile(outdir,'noise_by_ctrl.csv'));
save(fullfile(outdir,'noise_models.mat'), 'NF', 'files');
end

% ------------------------------------------------------------------------
function E = errs(file, ds, win)
% estimator errors at controller rate (as estimation_errors): e_j(k) = meas_j(k) - truth_j(k - ds(j))
D = load(file); R = D.R; x = double(D.Lx);
k = round(R.tc*1000); i = k + 1; sel = R.tc >= win(1) & R.tc < win(2);
k = k(sel); i = i(sel); meas = double(R.meas(sel,:));
phi = x(:,7); th = x(:,8); p = x(:,10); q = x(:,11); r = x(:,12);
ed = [p + sin(phi).*tan(th).*q + cos(phi).*tan(th).*r, cos(phi).*q - sin(phi).*r, (sin(phi).*q + cos(phi).*r)./cos(th)];
truth = [x(:,1:9), ed]; e = zeros(numel(i),12);
for j = 1:12, e(:,j) = meas(:,j) - truth(i - ds(j), j); end
E = struct('e', e, 'seg', cumsum([1; diff(k) ~= 1]), 'ctrl', D.c, 'seed', D.seed);
end

function Y = classdata(dd, cl)
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
Y = classdata(dd, cl); n = 0; s = 0; nd = 0; sd = 0;
for i = 1:numel(Y)
    y = Y{i} - mean(Y{i},1); d = diff(Y{i});
    s = s + sum(y.^2,1); n = n + size(y,1); sd = sd + sum(d.^2,1); nd = nd + size(d,1);
end
v = s/n; vd = sd/nd;
end

function [V, bic] = varfit(Y, n)
p = size(Y{1},2); G = zeros(n*p); H = zeros(n*p, p); N = 0;
for i = 1:numel(Y)
    y = Y{i}; m = size(y,1); Z = zeros(m-n, n*p);
    for l = 1:n, Z(:, (l-1)*p+(1:p)) = y(n+1-l:m-l, :); end
    G = G + Z'*Z; H = H + Z'*y(n+1:m,:); N = N + m - n;
end
B = G\H; Sig = 0;
for i = 1:numel(Y)
    y = Y{i}; m = size(y,1); Z = zeros(m-n, n*p);
    for l = 1:n, Z(:, (l-1)*p+(1:p)) = y(n+1-l:m-l, :); end
    Ee = y(n+1:m,:) - Z*B; Sig = Sig + Ee'*Ee;
end
Sig = Sig/N; V = struct('A', B', 'Sig', Sig, 'n', n, 'p', p);
bic = N*log(det(Sig)) + n*p^2*log(N);
end

function m = var2filter(V, sc)
p = V.p; n = V.n; Lm = chol(V.Sig, 'lower'); S = diag(sc);
Af = [V.A; eye(p*(n-1)), zeros(p*(n-1), p)]; Bf = [Lm; zeros(p*(n-1), p)];
m = struct('Af', Af, 'Bf', Bf, 'Cf', S*V.A, 'Df', S*Lm, 'order', NaN, 'ctrls', [], 'scale', sc);
end

function [v, vd] = filterstats(m)
P = dlyap(m.Af, m.Bf*m.Bf'); V = m.Cf*P*m.Cf' + m.Df*m.Df';
H = m.Cf*(m.Af - eye(size(m.Af,1))); Q = m.Cf*m.Bf - m.Df; Vd = H*P*H' + Q*Q' + m.Df*m.Df';
v = diag(V)'; vd = diag(Vd)';
end
