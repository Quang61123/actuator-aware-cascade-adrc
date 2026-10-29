function run_sitl_redesign(phase, worker, nworkers, efIds)
%RUN_SITL_REDESIGN Redesign PX4 SITL runs, executed exactly from the fixed plan
%   ../results/analysis/run_plan.csv (make_run_plan). Same SITL setup,
%   logging and metrics as run_sitl_estimator (EKF2 input, debug_array at 1 kHz).
%   phase: 'P1', 'P2', 'P3' or 'P4' (rows of phase 'P1b_unused' must not be run); efIds:
%   controller ids of the final -E designs ([17 18 19]) for the "EF1..EF3" rows of phase P2.
%   Each run is attempted at most twice only on an infrastructure error
%   (no ulog written), as in run_sitl_estimator; every run with a ulog is kept and
%   reported whatever its outcome. Output: ../results/analysis/runs/<tag>.mat
%   (M, X, run) and <work dir>\redesign_data\<tag>.mat (R, Lx).
if nargin < 2, worker = 0; nworkers = 1; end
if nargin < 4, efIds = [17 18 19]; end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
W = px4sitl_work_dir(); UW = fullfile(W,'ulog_redesign'); DW = fullfile(W,'redesign_data');
if ~exist(UW,'dir'), mkdir(UW); end
if ~exist(DW,'dir'), mkdir(DW); end
RA = fullfile(here,'..','results','analysis'); RR = fullfile(RA,'runs'); if ~exist(RR,'dir'), mkdir(RR); end
fp = fullfile(RA,'run_plan.csv'); if strcmp(phase,'P3'), fp = fullfile(RA,'run_plan_commanded.csv'); end   % P3: additional plan (make_run_plan_commanded)
if strcmp(phase,'P4'), fp = fullfile(RA,'run_plan_setpoint.csv'); end   % P4: PID-Budget (b*) (make_run_plan_setpoint)
io = detectImportOptions(fp);
io = setvartype(io, {'Phase','Group','CtrlId','Scenario','Split'}, 'string'); P = readtable(fp, io);
P.Split(ismissing(P.Split)) = "";
P = P(P.Phase == string(phase), :);
for r = 1:height(P)
    if mod(r-1, nworkers) ~= worker, continue; end
    cid = P.CtrlId(r); if startsWith(cid, "EF"), c = efIds(str2double(extractAfter(cid, 2))); else, c = str2double(cid); end
    hv = P.Scenario(r) == "hover"; scen = 6; if P.Scenario(r) == "S5", scen = 5; end
    seed = P.Seed(r); win = [P.WinStart_s(r) P.WinEnd_s(r)]; Tend = P.Tend_s(r);
    tag = sprintf('%s_%s_c%d_seed%d_T%d', P.Group(r), P.Scenario(r), c, seed, Tend);
    fout = fullfile(RR, [tag '.mat']); if exist(fout,'file'), continue; end
    opt = struct('seed', seed, 'tag', tag, 'instance', worker, 'keepUlog', true, 'hover', hv, 'raw', P.Group(r) == "PBR", 'Tend', Tend);
    ok = false;
    for attempt = 1:2
        try
            L = sitl_run(c, scen, opt); ok = ~isempty(L.ulog) && exist(L.ulog,'file'); if ok, break; end
        catch e
            logline(sprintf('[%s] attempt %d failed: %s', tag, attempt, e.message));
        end
    end
    if ~ok, logline(sprintf('[%s] RUN_FAILED no ulog', tag)); save(fout, 'tag'); continue; end
    uf = fullfile(UW,[tag '.ulg']); movefile(L.ulog, uf, 'f');
    try
        R = sitl_log_extract(uf);
    catch e
        logline(sprintf('[%s] extract failed: %s', tag, e.message)); R = [];
    end
    pv = quad_params(); X = struct();
    if ~isempty(R)
        kw = R.tc >= win(1) & R.tc < win(2); okw = kw & R.ok; okw(1) = false; ip = find(okw);
        dW = pv(11)*double(R.wn(ip,:) - R.wn(ip-1,:)); dU = double(R.ua(ip,:) - R.ua(ip-1,:));
        X = struct('J_unq', sqrt(mean(dW(:).^2)), 'sd_du', std(dU,1,1), 'sd_u', std(double(R.ua(kw,:)),1,1), ...
            'n_pairs', numel(ip), 'n_expected', round(diff(win)*1000));
        Lx = single(L.x); Lt = L.t; tArm = L.tArm; %#ok<NASGU>
        save(fullfile(DW,[tag '.mat']), 'R', 'Lx', 'Lt', 'tArm', 'c', 'seed', 'hv', 'win', '-v7.3');
    end
    kq = L.t >= win(1) & L.t < win(2); dwq = diff(L.wcmd(kq,:)); X.J_q = sqrt(mean(dwq(:).^2));
    M = run_metrics(L, win, pv);
    run = struct('plan_id', P.RunID(r), 'group', P.Group(r), 'ctrl', c, 'hover', hv, 'scen', scen, 'seed', seed, 'Tend', Tend, ...
        'split', P.Split(r), 'tag', tag, 'win', win, 'miss', L.miss, 'wall', L.wall, 'diverged', L.diverged, 'info', L.info); %#ok<NASGU>
    save(fout, 'M', 'X', 'run');
    logline(sprintf('[redesign %s] rmse %.1f mm, J_unq %.4f, miss %d, div %d', tag, M.rmse3*1e3, getfield_or(X,'J_unq'), L.miss, L.diverged));
end
logline(sprintf('RUN_DONE phase %s worker %d', phase, worker));
end

function v = getfield_or(X, f), if isfield(X, f), v = X.(f); else, v = NaN; end, end

function logline(s)
fid = fopen([fullfile(px4sitl_work_dir(),'w') getenv('SITL_WORKER') '.log'],'a');
fprintf(fid,'%s %s\n', char(datetime('now','Format','HH:mm:ss')), s); fclose(fid);
end
