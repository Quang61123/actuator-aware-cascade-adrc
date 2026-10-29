function L = sitl_run(ctrlId, scen, opt)
%SITL_RUN One PX4 software-in-the-loop run of the controller.
%   The nonlinear 16-state plant (quad_dynamics, RK4, 1 ms) runs
%   here in lockstep with PX4 v1.16.1 SITL (WSL). Every 1 ms the plant sends
%   HIL_SENSOR (IMU; mag/baro every 10 ms) and every 10 ms an ODOMETRY
%   message (motion-capture-like position, velocity, attitude, 5 ms latency).
%   PX4 runs EKF2 and the code-generated controller (px4_simulink_app) and
%   returns HIL_ACTUATOR_CONTROLS; motor i command = u_i * w_max.
%   ctrlId: 1 PID-Budget, 2 AA-Budget, 3 RO-HalfGain, 4 AA-HalfGain,
%           5 A1-Commanded, 6 A2-Unaligned, 7 A3-NoIncrement,
%           8 PID-Budget-S, 9 AA-Budget-S, 10 AA-HalfGain-S; further ids in create_sitl_controller.
%   scen: 6 set-point with turbulence (40 s), 5 lemniscate (45 s).
%   opt (optional): seed, Tend, warm [s], tag, pvt (true plant), verbose,
%   instance, keepUlog, hover, raw (diagnostic
%   estimator bypass; every 1 ms an additional MOCAP ODOMETRY message carries
%   the sensor model (sensor_model), read by the controller instead of EKF2).
%   L: t, x (Nx16), r (Nx10), wcmd (Nx4), miss, wall, armed time, info.
if nargin < 3, opt = struct(); end
here = fileparts(mfilename('fullpath'));
simRoot = fullfile(here,'..','MATLAB_simulation'); addpath(here, simRoot);
[W, Wwsl] = px4sitl_work_dir();
cfg = sim_config(); Gs = load(fullfile(cfg.results,'design_gains.mat'),'S'); S = Gs.S;
Ts = 1e-3;
seed = getf(opt,'seed', cfg.nominalSeed + 1);
hover = getf(opt,'hover', false);         % noise-only hover (scen = 0)
if hover, scen = 0; end
Tend = getf(opt,'Tend', 40 + 5*(scen == 5) + 20*hover);
warm = getf(opt,'warm', 12);              % pre-arm EKF2 warm-up [s, simulated]
tag  = getf(opt,'tag', sprintf('c%d_s%d_%d', ctrlId, scen, seed));
verbose = getf(opt,'verbose', true);
inst = getf(opt,'instance', 0);           % PX4 instance (parallel campaigns)
raw  = double(getf(opt,'raw', false));    % diagnostic estimator bypass (controller reads the raw 1-kHz measurement)
pvt = quad_params();
pvt = getf(opt,'pvt', pvt);
N = round(Tend/Ts);
if hover   % as numerical_checks: make_exo(2,...), no bias, no wind
    E = make_exo(2, seed, Tend, Ts, S); E.bias(:) = 0; E.wind(:) = 0; noisy = true;
else
    E = make_exo(scen, seed, Tend, Ts, S); noisy = scen >= 2;
end
rs = RandStream('mt19937ar','Seed',seed + 100000);
sigA = 0.05*noisy; sigM = 0.002*noisy; sigB = 0.01*noisy;   % accel [m/s^2], mag [G], baro [hPa]
NA = randn(rs, 7, N);
Bned = [0.40; 0; 0.10];                    % assumed earth field [G], zero declination
evLag = 5;                                 % EV latency [samples] (EKF2_EV_DELAY = 5 ms)

% ---------- start relay + PX4 in WSL ----------
prepare_wsl(here, W);
wsl(sprintf(['bash ' Wwsl '/px4_sitl.sh start %d %d %s %d %d'], ctrlId, scen, tag, inst, raw), false);
pause(1);
c = []; t0 = tic;   % .NET socket (NoDelay, low latency on Windows)
while isempty(c) && toc(t0) < 20
    try, c = System.Net.Sockets.TcpClient(); c.NoDelay = true; c.Connect('127.0.0.1', 4600 + inst); catch, c = []; pause(0.5); end
end
assert(~isempty(c), 'sitl_run: cannot connect to relay');
ns = c.GetStream(); rb = NET.createArray('System.Byte', 4096);
cleanup = onCleanup(@() stop_all(c, inst));

% ---------- logs ----------
L.t = (0:N-1)'*Ts; L.x = nan(N,16); L.r = nan(N,10); L.wcmd = nan(N,4);
L.aref = nan(N,3); L.miss = 0; L.ctrlId = ctrlId; L.scen = scen; L.seed = seed;
x = zeros(16,1); x(13:16) = sqrt(pvt(1)*pvt(2)/(4*pvt(4)));
xh = repmat(x,1,evLag+1); nh = zeros(12,evLag+1); bh = zeros(3,evLag+1);   % EV history (col 1 newest)
wc = x(13:16); buf = zeros(1,0,'uint8'); seq = 0;
t_us = 0; step = 0; phase = 0; k = 0; armSent = false; tArm = NaN; lastU = zeros(1,16);
tw = tic; L.diverged = false; lastT = 0; acts = [];
while true
    step = step + 1; t_us = t_us + 1000;
    if phase == 1, i = k + 1; else, i = 1; end
    if phase == 1, n12 = E.noise(:,i); bias = E.bias(:,i); else, n12 = zeros(12,1); bias = zeros(3,1); end
    % ----- IMU (FRD) -----
    [Fd, Td] = disturbance_gen((i-1)*Ts, scen);
    if phase == 1, Fd = Fd + E.wind(:,i); xd = quad_dynamics(x, wc, Fd, Td, pvt); acc = xd(4:6);
    else, acc = zeros(3,1); end
    [~,~,~,gyr,Rn] = state_to_px4(x);
    Rflu = Rn; Rflu(:,2:3) = -Rflu(:,2:3); Rflu(2:3,:) = -Rflu(2:3,:);    % FLU->NWU
    f = Rflu'*(acc + [0;0;pvt(2)]);
    accF = [f(1); -f(2); -f(3)] + sigA*NA(1:3,i);
    gyrF = gyr + S.sig(4)*[n12(10); -n12(11); -n12(12)] + [bias(1); -bias(2); -bias(3)];
    full = mod(step,10) == 0;
    magF = Rn'*Bned + sigM*NA(4:6,i);
    alt = 10 + x(3); pabs = 1013.25*(1 - 2.25577e-5*alt)^5.25588 + sigB*NA(7,i);
    fields = 63; if full, fields = 63 + 448 + 512 + 2048 + 4096; end
    pl = [typecast(uint64(t_us),'uint8'), typecast(single([accF; gyrF; magF; pabs; 0; alt; 20]'),'uint8'), ...
          typecast(uint32(fields),'uint8'), uint8(0)];
    out = mav_frame(107, pl, 108, seq); seq = seq + 1;
    % ----- external vision (every 10 ms, 5 ms latency) -----
    xh = [x, xh(:,1:end-1)]; nh = [n12, nh(:,1:end-1)]; bh = [bias, bh(:,1:end-1)];
    if raw   % diagnostic estimator bypass: output of sensor_model every 1 ms (n_d = 5)
        xr = xh(:,end); nr = nh(:,end); br = bh(:,end);
        xr(1:3) = xr(1:3) + S.sig(1)*nr(1:3); xr(4:6) = xr(4:6) + S.sig(2)*nr(4:6);
        xr(7:9) = xr(7:9) + S.sig(3)*nr(7:9); xr(10:12) = xr(10:12) + br + S.sig(4)*nr(10:12);
        [pr, vr, qr, gr] = state_to_px4(xr);
        pl = [typecast(uint64(t_us - 1000*evLag),'uint8'), ...
              typecast(single([pr; qr; vr; gr; zeros(42,1)]'),'uint8'), uint8([1 1 0 6 100])];
        out = [mav_frame(331, pl, 91, seq), out]; seq = seq + 1;   % MOCAP -> vehicle_mocap_odometry, sent before HIL_SENSOR
    end
    if full
        xe = xh(:,end); ne = nh(:,end);
        xe(1:3) = xe(1:3) + S.sig(1)*ne(1:3); xe(4:6) = xe(4:6) + S.sig(2)*ne(4:6);
        xe(7:9) = xe(7:9) + S.sig(3)*ne(7:9);
        [pe, ve, qe] = state_to_px4(xe);
        pc = zeros(1,21); pc([1 7 12]) = max(S.sig(1),1e-3)^2; pc([16 19 21]) = max(S.sig(3),1e-3)^2;
        vc = zeros(1,21); vc([1 7 12]) = max(S.sig(2),1e-3)^2;
        pl = [typecast(uint64(t_us - 1000*evLag),'uint8'), ...
              typecast(single([pe; qe; ve; nan(3,1); pc'; vc']'),'uint8'), uint8([1 1 0 2 100])];
        out = [out, mav_frame(331, pl, 91, seq)]; seq = seq + 1; %#ok<AGROW>
    end
    ns.Write(out, 0, numel(out));
    % ----- wait for the lockstep actuator response -----
    if phase == 0, tmo = 0.0015; else, tmo = 0.5; end
    got = false; ta = tic;
    while true
        nb = c.Available;
        if nb > 0
            nr = ns.Read(rb, 0, min(nb,4096)); bb = uint8(rb); buf = [buf, bb(1:nr)]; %#ok<AGROW>
            [acts, buf] = mav_parse_act(buf);
            for a = 1:numel(acts)
                lastU = acts(a).u; lastMode = acts(a).mode; lastT = acts(a).t;
                if acts(a).t >= t_us, got = true; end
            end
        end
        if got || toc(ta) > tmo, break; end
    end
    if phase == 0
        if step == round(warm/Ts) && ~armSent
            wsl(sprintf(['bash ' Wwsl '/px4_sitl.sh arm %d'], inst), true); armSent = true;
            if verbose, logmsg(W, '[%s] arm command at %.1f s (wall %.1f s)\n', tag, step*Ts, toc(tw)); end
        end
        if got && exist('lastMode','var') && bitand(lastMode,128) && any(lastU(1:4) > 0)
            phase = 1; tArm = step*Ts;
            if verbose, logmsg(W, '[%s] armed, flight starts (PX4 time %.3f s)\n', tag, tArm); end
        elseif step*Ts > warm + 30
            error('sitl_run:arm','[%s] vehicle did not arm', tag);
        end
        if phase == 0, continue; end
    end
    % ----- flight: apply command, log, integrate -----
    if ~got, L.miss = L.miss + 1; end
    if bitand(lastMode,128), wc = lastU(1:4)'*pvt(11); else, wc = zeros(4,1); end
    wc = min(max(wc,0), pvt(11));
    L.x(i,:) = x'; if scen == 0, L.r(i,:) = 0; else, L.r(i,:) = reference_gen((i-1)*Ts, scen)'; end
    L.wcmd(i,:) = wc';
    [Fd, Td] = disturbance_gen((i-1)*Ts, scen); Fd = Fd + E.wind(:,i);
    h = Ts;
    k1 = quad_dynamics(x, wc, Fd, Td, pvt); k2 = quad_dynamics(x+h/2*k1, wc, Fd, Td, pvt);
    k3 = quad_dynamics(x+h/2*k2, wc, Fd, Td, pvt); k4 = quad_dynamics(x+h*k3, wc, Fd, Td, pvt);
    x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
    k = k + 1;
    if any(~isfinite(x)) || max(abs(x(1:3))) > 50
        L.diverged = true; if verbose, logmsg(W, '[%s] diverged at %.2f s\n', tag, k*Ts); end
        break
    end
    if k >= N, break; end
    if verbose && mod(k,10000) == 0, logmsg(W, '[%s] t = %5.1f s  wall %.0f s  miss %d\n', tag, k*Ts, toc(tw), L.miss); end
end
L.wall = toc(tw); L.tArm = tArm; L.tag = tag;
L.ulog = '';
if getf(opt,'keepUlog', false)
    [~,o] = system(sprintf(['wsl -d PX4-Ubuntu-22.04 -- bash ' Wwsl '/px4_sitl.sh lastlog %d'], inst));
    p = strtrim(o); j = strfind(p,'/home/'); 
    if ~isempty(j), p = p(j(1):end); L.ulog = fullfile(W, [tag '.ulg']); stop_all(c, inst); pause(0.5);
        copyfile(['\\wsl$\PX4-Ubuntu-22.04' strrep(p,'/','\')], L.ulog); end
end

L.info = struct('plant','16-state RK4 1 ms','px4','v1.16.1 SITL lockstep','ev','10 ms, 5 ms latency', ...
    'sigA',sigA,'sigM',sigM,'sigB',sigB,'evLag',evLag,'raw',raw);
end

% ------------------------------------------------------------------------
function logmsg(W, varargin)
fid = fopen(fullfile(W,sprintf('sitl_progress_%s.log', getenv('SITL_WORKER'))),'a'); fprintf(fid, varargin{:}); fclose(fid);
end

function v = getf(s, f, d)
if isfield(s,f), v = s.(f); else, v = d; end
end

function wsl(cmd, async)
f = fullfile(W, '_run.cmd');
if async
    system(['start "" /b wsl -d PX4-Ubuntu-22.04 -- ' cmd ' > NUL 2>&1']);
else
    system(['wsl -d PX4-Ubuntu-22.04 -- ' cmd]); %#ok<*NASGU>
end
end

function stop_all(c, inst)
[~, Wwsl] = px4sitl_work_dir();
try, c.Close(); catch, end
system(sprintf(['wsl -d PX4-Ubuntu-22.04 -- bash ' Wwsl '/px4_sitl.sh stop %d > NUL 2>&1'], inst));
end

function prepare_wsl(here, W)
[~, Wwsl] = px4sitl_work_dir();
persistent done
if ~isempty(done) || ~isempty(getenv('SITL_WORKER')), return; end   % parallel workers: files prepared beforehand
copyfile(fullfile(here,'wsl','*'), W, 'f');
system(['wsl -d PX4-Ubuntu-22.04 -- bash -c "sed -i ''s/\r$//'' ' Wwsl '/px4_sitl.sh ' Wwsl '/10099_ladrc_sitl ' Wwsl '/sim_relay.py"']);
done = true;
end
