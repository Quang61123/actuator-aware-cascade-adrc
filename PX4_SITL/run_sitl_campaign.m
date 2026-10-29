function T = run_sitl_campaign(parts, worker, nworkers)
%RUN_SITL_CAMPAIGN Main PX4 SITL campaign (resumable).
%   Controller ids: 1 PID-Budget, 2 AA-Budget, 3 RO-HalfGain, 4 AA-HalfGain,
%   5 A1-Commanded, 6 A2-Unaligned, 7 A3-NoIncrement, 8 PID-Budget-S,
%   9 AA-Budget-S, 10 AA-HalfGain-S. Gains: results/design_gains.mat.
%   parts: any of
%     'A' S6, primary designs 1-4, paired seeds 91001-91030 (as E1)
%     'B' S6, strict designs 8-10, seeds 91001-91030 (as E0-strict)
%     'C' S5 lemniscate, designs 1-4, seeds 91001-91010
%     'D' S6, ablations 5-7, seeds 91001-91010
%     'E' S6, plant rotor lag 75/90/105/120 ms (n_d = 5), designs 1, 2, 4,
%         seeds 92001-92003 (internal model fixed at 30 ms, as E7)
%     'F' repeatability: AA-Budget, S6, seed 91001, two more runs
%   worker/nworkers (optional): this process runs plan rows with
%   mod(row-1, nworkers) == worker on PX4 instance 'worker' (parallel runs).
%   Each run: results/runs/<tag>.mat (metrics + 10-ms trajectories);
%   full-rate logs in <work dir>\full_main.
here = fileparts(mfilename('fullpath'));
simRoot = fullfile(here,'..','MATLAB_simulation'); addpath(here, simRoot);
W = px4sitl_work_dir(); FW = fullfile(W,'full_main'); if ~exist(FW,'dir'), mkdir(FW); end
R = fullfile(here,'results','runs'); if ~exist(R,'dir'), mkdir(R); end
plan = zeros(0,5);   % [ctrl scen seed tau_plant_ms rep]
if contains(parts,'A'), for q = 1:30, for c = 1:4,  plan(end+1,:) = [c 6 91000+q 0 0]; end, end, end %#ok<AGROW>
if contains(parts,'B'), for q = 1:30, for c = 8:10, plan(end+1,:) = [c 6 91000+q 0 0]; end, end, end %#ok<AGROW>
if contains(parts,'C'), for q = 1:10, for c = 1:4,  plan(end+1,:) = [c 5 91000+q 0 0]; end, end, end %#ok<AGROW>
if contains(parts,'D'), for q = 1:10, for c = 5:7,  plan(end+1,:) = [c 6 91000+q 0 0]; end, end, end %#ok<AGROW>
if contains(parts,'E')
    for tm = [75 90 105 120], for q = 1:3, for c = [1 2 4], plan(end+1,:) = [c 6 92000+q tm 0]; end, end, end %#ok<AGROW>
end
if contains(parts,'F'), for rep = 1:2, plan(end+1,:) = [2 6 91001 0 rep]; end, end %#ok<AGROW>
if nargin < 2, worker = 0; nworkers = 1; end
for r = 1:size(plan,1)
    if mod(r-1, nworkers) ~= worker, continue; end
    c = plan(r,1); s = plan(r,2); seed = plan(r,3); tm = plan(r,4); rep = plan(r,5);
    tag = sprintf('c%d_s%d_seed%d_tau%d_rep%d', c, s, seed, tm, rep);
    fout = fullfile(R, [tag '.mat']);
    if exist(fout,'file'), continue; end
    opt = struct('seed', seed, 'tag', tag, 'instance', worker, 'keepUlog', (c == 2 && s == 6 && seed == 91002 && rep == 0 && tm == 0));
    if tm > 0, opt.pvt = quad_params([], tm/1000, 5); end
    ok = false;
    for attempt = 1:2
        try
            L = sitl_run(c, s, opt); ok = true; break
        catch e
            logline(W, sprintf('[%s] attempt %d failed: %s', tag, attempt, e.message));
        end
    end
    if ~ok, continue; end
    save(fullfile(FW,[tag '.mat']), 'L', '-v7.3');
    M = sitl_metrics(L);
    D.t = L.t(1:10:end); D.x = single(L.x(1:10:end,:)); D.r = single(L.r(1:10:end,:));
    D.wcmd = single(L.wcmd(1:10:end,:));
    run = struct('ctrl', c, 'scen', s, 'seed', seed, 'tau_ms', tm, 'rep', rep, 'tag', tag, ...
        'miss', L.miss, 'wall', L.wall, 'diverged', L.diverged, 'ulog', L.ulog, 'info', L.info); %#ok<NASGU>
    save(fout, 'M', 'D', 'run');
    logline(W, sprintf('[main %s] done: rmse %.1f mm, inc %.2f, div %d, miss %d, wall %.0f s', tag, M.rmse3*1e3, M.chat_cmd, L.diverged, L.miss, L.wall));
end
T = sitl_collect();
end

function logline(W, s)
fid = fopen(fullfile(W,['campaign_main_w' getenv('SITL_WORKER') '.log']),'a'); fprintf(fid,'%s %s\n', char(datetime('now','Format','HH:mm:ss')), s); fclose(fid);
end
