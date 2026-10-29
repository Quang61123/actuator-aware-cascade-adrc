function run_sitl_estimator(parts, worker, nworkers)
%RUN_SITL_ESTIMATOR PX4 SITL runs with controller-rate logging (estimator stack), resumable.
%   Same SITL configuration as run_sitl_campaign (PX4 v1.16.1 lockstep, 1-kHz IMU,
%   100-Hz external vision with 5-ms latency, nominal sensor noise), gains and
%   control logic unchanged; the controller additionally logs, at every
%   1-ms step, the measurement vector it uses (debug_array), and the plant
%   state is kept, so that estimator errors can be formed at controller rate.
%   'H' noise-only hover, EKF2 input, 6 designs, seeds 96001-96010, 60 s
%   'S' S6, EKF2 input, PID-Budget/AA-Budget/AA-HalfGain, seeds 91001-91003
%       (effective estimator delay)
%   'C' calibration of the controller/plant time alignment: raw-measurement
%       input (estimator bypass), AA-Budget, hover, seed 96001
%   Output: ../results/estimator_noise/runs/<tag>.mat (M, X, run) and the
%   controller-rate data <work dir>\est_data\<tag>.mat (R, Lx).
if nargin < 2, worker = 0; nworkers = 1; end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
W = px4sitl_work_dir(); UW = fullfile(W,'ulog_est'); DW = fullfile(W,'est_data');
if ~exist(UW,'dir'), mkdir(UW); end
if ~exist(DW,'dir'), mkdir(DW); end
RR = fullfile(here,'..','results','estimator_noise','runs'); if ~exist(RR,'dir'), mkdir(RR); end
%   Equal-activity designs for the estimator-fed loop ( ids 11 AA-Budget-E,
%   12 AA-HalfGain-E, see design_estimator_initial):
%   'E' S6, PID-Budget and the two -E designs, seeds 91001-91030
%   'F' noise-only hover, the two -E designs, seeds 96001-96010
%   'G' S5 (lemniscate), PID-Budget and the two -E designs, seeds 91001-91010
plan = zeros(0,6);   % [ctrl hover seed raw scen group(0 base, 1 -E designs, 2 -E2 designs)]
if contains(parts,'C'), plan(end+1,:) = [2 1 96001 1 6 0]; end
if contains(parts,'H'), for q = 1:10, for c = [1 2 4 8 9 10], plan(end+1,:) = [c 1 96000+q 0 6 0]; end, end, end %#ok<AGROW>
if contains(parts,'S'), for q = 1:3, for c = [1 2 4], plan(end+1,:) = [c 0 91000+q 0 6 0]; end, end, end %#ok<AGROW>
if contains(parts,'E'), for q = 1:30, for c = [1 11 12], plan(end+1,:) = [c 0 91000+q 0 6 1]; end, end, end %#ok<AGROW>
if contains(parts,'F'), for q = 1:10, for c = [11 12], plan(end+1,:) = [c 1 96000+q 0 6 1]; end, end, end %#ok<AGROW>
if contains(parts,'G'), for q = 1:10, for c = [1 11 12], plan(end+1,:) = [c 0 91000+q 0 5 1]; end, end, end %#ok<AGROW>
%   Yaw-constrained designs (ids 13 AA-Budget-E2, 14 AA-HalfGain-E2,
%   results/yaw_ceiling/yaw_ceiling_designs.mat; output to ../results/yaw_ceiling/runs):
%   'K' S6, the two -E2 designs, seeds 91001-91030
%   'L' noise-only hover, the two -E2 designs, seeds 96001-96010
if contains(parts,'K'), for q = 1:30, for c = [13 14], plan(end+1,:) = [c 0 91000+q 0 6 2]; end, end, end %#ok<AGROW>
if contains(parts,'L'), for q = 1:10, for c = [13 14], plan(end+1,:) = [c 1 96000+q 0 6 2]; end, end, end %#ok<AGROW>
%   'M' S5 (lemniscate), the two -E2 designs, seeds 91001-91010
if contains(parts,'M'), for q = 1:10, for c = [13 14], plan(end+1,:) = [c 0 91000+q 0 5 2]; end, end, end %#ok<AGROW>
RRyaw = fullfile(here,'..','results','yaw_ceiling','runs'); if ~exist(RRyaw,'dir'), mkdir(RRyaw); end
RRest = RR;
for r = 1:size(plan,1)
    if mod(r-1, nworkers) ~= worker, continue; end
    c = plan(r,1); hv = plan(r,2); seed = plan(r,3); raw = plan(r,4); scen = plan(r,5);
    pre = 'est'; if raw, pre = 'est_cal'; end
    if plan(r,6) == 1, pre = 'est_E'; end
    RR = RRest; if plan(r,6) == 2, pre = 'yaw_ceiling'; RR = RRyaw; end
    if hv, tag = sprintf('%s_hover_c%d_seed%d', pre, c, seed); win = [20 60];
    elseif scen == 5, tag = sprintf('%s_s5_c%d_seed%d', pre, c, seed); win = [5 45];
    else,  tag = sprintf('%s_s6_c%d_seed%d', pre, c, seed);    win = [10 40]; end
    fout = fullfile(RR, [tag '.mat']); if exist(fout,'file'), continue; end
    opt = struct('seed', seed, 'tag', tag, 'instance', worker, 'keepUlog', true, 'hover', logical(hv), 'raw', logical(raw));
    ok = false;
    for attempt = 1:2
        try
            L = sitl_run(c, scen, opt); ok = ~isempty(L.ulog) && exist(L.ulog,'file'); if ok, break; end
        catch e
            logline(sprintf('[%s] attempt %d failed: %s', tag, attempt, e.message));
        end
    end
    if ~ok, logline(sprintf('[%s] no ulog', tag)); continue; end
    uf = fullfile(UW,[tag '.ulg']); movefile(L.ulog, uf, 'f');
    try
        R = sitl_log_extract(uf);
    catch e
        logline(sprintf('[%s] extract failed: %s', tag, e.message)); continue
    end
    % activity in the metric window from the controller-rate log (consecutive steps only)
    pv = quad_params(); kw = R.tc >= win(1) & R.tc < win(2); okw = kw & R.ok; okw(1) = false;
    ip = find(okw); dW = pv(11)*double(R.wn(ip,:) - R.wn(ip-1,:)); dU = double(R.ua(ip,:) - R.ua(ip-1,:));
    X = struct('J_unq', sqrt(mean(dW(:).^2)), 'sd_du', std(dU,1,1), 'sd_u', std(double(R.ua(kw,:)),1,1), ...
        'n_pairs', numel(ip), 'n_expected', round(diff(win)*1000));
    kq = L.t >= win(1) & L.t < win(2); dwq = diff(L.wcmd(kq,:)); X.J_q = sqrt(mean(dwq(:).^2));
    Lx = single(L.x); Lt = L.t; tArm = L.tArm; %#ok<NASGU>
    save(fullfile(DW,[tag '.mat']), 'R', 'Lx', 'Lt', 'tArm', 'c', 'seed', 'hv', 'raw', 'win', '-v7.3');
    M = run_metrics(L, win, pv);
    run = struct('ctrl', c, 'hover', hv, 'raw', raw, 'scen', scen, 'seed', seed, 'tag', tag, 'win', win, ...
        'miss', L.miss, 'wall', L.wall, 'diverged', L.diverged, 'info', L.info, 'n_dbg', R.n, 'frac_dbg_ok', R.frac_ok); %#ok<NASGU>
    save(fout, 'M', 'X', 'run');
    logline(sprintf('[est %s] rmse %.1f mm, J_unq %.3f, sd_du %.4f/%.4f, dbg %d (ok %.3f), miss %d', ...
        tag, M.rmse3*1e3, X.J_unq, X.sd_du(1), X.sd_du(2), R.n, R.frac_ok, L.miss));
end
end

function logline(s)
fid = fopen([fullfile(px4sitl_work_dir(),'est_w') getenv('SITL_WORKER') '.log'],'a');
fprintf(fid,'%s %s\n', char(datetime('now','Format','HH:mm:ss')), s); fclose(fid);
end
