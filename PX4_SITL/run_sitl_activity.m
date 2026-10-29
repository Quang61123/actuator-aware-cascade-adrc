function run_sitl_activity(parts, worker, nworkers)
%RUN_SITL_ACTIVITY PX4 SITL activity-check runs, resumable.
%   'S' S6, PID-Budget, AA-Budget, AA-HalfGain and the three -S designs,
%       paired seeds 91001-91010 (window 10-40 s)
%   'H' noise-only hover at the origin (no bias, wind or disturbance), same
%       designs, seeds 96001-96010, 60 s (window 20-60 s), as in the MATLAB hover check (hover_per_seed.csv)
%   'R','Q' as 'S','H' with the estimator bypass: the controller reads the
%       sensor model of sensor_model (nominal white noise, gyro bias,
%       5-ms delay, 1 kHz) instead of the EKF2 estimates (diagnostic)
%   Gains, internal models and control logic are unchanged; the
%   controller only exports the unquantized command and u_phi, u_theta.
%   Output: ../results/sitl_activity/runs/<tag>.mat (M, X, run); ulogs in
%   <work dir>\ulog_act (local, large).
if nargin < 2, worker = 0; nworkers = 1; end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
W = px4sitl_work_dir(); UW = fullfile(W,'ulog_act'); if ~exist(UW,'dir'), mkdir(UW); end
R = fullfile(here,'..','results','sitl_activity','runs'); if ~exist(R,'dir'), mkdir(R); end
ids = [1 2 4 8 9 10];
plan = zeros(0,4);   % [ctrl hover seed raw]
for pr = {'S',0,0; 'H',1,0; 'R',0,1; 'Q',1,1}'
    if ~contains(parts, pr{1}), continue; end
    for q = 1:10, for c = ids, plan(end+1,:) = [c pr{2} 91000+5000*pr{2}+q pr{3}]; end, end %#ok<AGROW>
end
for r = 1:size(plan,1)
    if mod(r-1, nworkers) ~= worker, continue; end
    c = plan(r,1); hv = plan(r,2); seed = plan(r,3); raw = plan(r,4);
    pre = 'act'; if raw, pre = 'act_raw'; end
    if hv, tag = sprintf('%s_hover_c%d_seed%d', pre, c, seed); win = [20 60];
    else,  tag = sprintf('%s_s6_c%d_seed%d', pre, c, seed);    win = [10 40]; end
    fout = fullfile(R, [tag '.mat']); if exist(fout,'file'), continue; end
    opt = struct('seed', seed, 'tag', tag, 'instance', worker, 'keepUlog', true, 'hover', logical(hv), 'raw', logical(raw));
    ok = false;
    for attempt = 1:2
        try
            L = sitl_run(c, 6, opt); ok = ~isempty(L.ulog) && exist(L.ulog,'file'); if ok, break; end
        catch e
            logline(sprintf('[%s] attempt %d failed: %s', tag, attempt, e.message));
        end
    end
    if ~ok, logline(sprintf('[%s] no ulog', tag)); continue; end
    movefile(L.ulog, fullfile(UW,[tag '.ulg']), 'f');
    try
        X = sitl_activity_extract(L, fullfile(UW,[tag '.ulg']), win);
    catch e
        logline(sprintf('[%s] extract failed: %s', tag, e.message)); continue
    end
    pv = quad_params(); M = run_metrics(L, win, pv);
    run = struct('ctrl', c, 'hover', hv, 'raw', raw, 'seed', seed, 'tag', tag, 'win', win, ...
        'miss', L.miss, 'wall', L.wall, 'diverged', L.diverged, 'info', L.info); %#ok<NASGU>
    save(fout, 'M', 'X', 'run');
    logline(sprintf('[act %s] rmse %.1f mm, J_q %.3f, J_unq %.3f, sd_du roll %.3f pitch %.3f, motors %d/%d, miss %d', ...
        tag, M.rmse3*1e3, X.J_q, X.J_unq, X.sd_du(1), X.sd_du(2), X.n_motor_samples, X.n_expected, L.miss));
end
end

function logline(s)
fid = fopen([fullfile(px4sitl_work_dir(),'act_w') getenv('SITL_WORKER') '.log'],'a');
fprintf(fid,'%s %s\n', char(datetime('now','Format','HH:mm:ss')), s); fclose(fid);
end
