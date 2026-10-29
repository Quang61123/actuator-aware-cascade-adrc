function X = sitl_activity_extract(L, ulogFile, win)
%SITL_ACTIVITY_EXTRACT Activity signals from the PX4 ulog of one SITL run.
%   actuator_motors (logged at every update, SDLOG_PROFILE high-rate bit):
%     control(1:4)*w_max = rotor-speed command BEFORE output quantization,
%     control(5:8)       = limited virtual controls u_phi, u_theta, u_psi, u_z,
%   vehicle_local_position / vehicle_attitude = EKF2 estimates, compared with
%   the plant state L.x. win = metric window [s] after arming.
%   Returns rotor-command increments J (on unquantized and quantized
%   commands), virtual-control statistics and estimator-error statistics.
pv = quad_params(); wmax = pv(11);
u = ulogreader(ulogFile);
d = readTopicMsgs(u,'TopicNames',{'actuator_motors','vehicle_local_position','vehicle_attitude'});
nm = string(d.TopicNames);
A = d.TopicMessages{nm == "actuator_motors"};
P = d.TopicMessages{nm == "vehicle_local_position"};
Q = d.TopicMessages{nm == "vehicle_attitude"};
% ---- actuator_motors in the window ----
ta = seconds(A.timestamp) - L.tArm; c = double(A.control);
ok = ta >= win(1) & ta < win(2) & c(:,10) > 0;            % c(10) = controller id (armed samples)
ta = ta(ok); c = c(ok,:);
[ta, iu] = unique(ta); c = c(iu,:);
X.n_motor_samples = numel(ta); X.n_expected = round(diff(win)*1000);
X.max_gap_ms = 1e3*max(diff(ta));
cons = round(diff(c(:,9))*1000) == 1;                     % consecutive controller steps only (logger may drop samples)
X.frac_consecutive = mean(cons);
wu = wmax*c(:,1:4); dwu = diff(wu); dwu = dwu(cons,:);
X.J_unq = sqrt(mean(dwu(:).^2));                            % rotor-command increment J, unquantized command
kq = L.t >= win(1) & L.t < win(2); dwq = diff(L.wcmd(kq,:));
X.J_q = sqrt(mean(dwq(:).^2));                              % same formula, quantized command received by the plant
X.frac_q_nonzero = mean(dwq(:) ~= 0);
V = c(:,[5 6 7 8]);                                         % [u_phi u_theta u_psi u_z] (= ua(4),ua(5),ua(6),ua(3) of ctrl_law)
X.sd_u = std(V,1,1);
dV = diff(V); dV = dV(cons,:);
X.sd_du = std(dV,1,1);                                 % increment std, same units as the MATLAB hover check (sigma_Delta u)
X.var_u = X.sd_u.^2; X.var_du = X.sd_du.^2;
X.ctrl_time_span = [min(c(:,9)) max(c(:,9))];
% ---- measurement as used by the controller (u(11) = x, u(12) = roll) ----
% compared with the plant state at controller time minus a lag of 0-15 ms
% (controller time k*Ts and plant time both start at the first armed step)
tc = c(:,9); lg = 0:15; ex = zeros(size(lg)); ef = ex;
for j = 1:numel(lg)
    xs = interp1(L.t, L.x(:,[1 7]), tc - lg(j)*1e-3);
    ex(j) = mean((c(:,11) - xs(:,1)).^2, 'omitnan'); ef(j) = mean((c(:,12) - xs(:,2)).^2, 'omitnan');
end
[~,j1] = min(ex); [~,j2] = min(ef); X.meas_lag_ms = [lg(j1) lg(j2)];
xs = interp1(L.t, L.x(:,[1 7]), tc - lg(j1)*1e-3); X.meas_err_sd_x = std(c(:,11) - xs(:,1), 1, 'omitnan');
xs = interp1(L.t, L.x(:,[1 7]), tc - lg(j2)*1e-3); X.meas_err_sd_phi = std(c(:,12) - xs(:,2), 1, 'omitnan');
dm = diff(c(:,[11 12])); X.meas_inc_sd = std(dm(cons,:), 1, 1);   % 1-ms increments of the measured x and roll
% ---- EKF2 estimation error (NED -> model frame) ----
tL = L.t(end);                                            % last plant sample (interpolation range)
tp = seconds(P.timestamp) - L.tArm; okp = tp >= win(1) & tp < min(win(2), tL); tp = tp(okp);
xt = interp1(L.t, L.x, tp);
ep = [P.x(okp) - xt(:,1), -P.y(okp) - xt(:,2), -P.z(okp) - xt(:,3)];
ev = [P.vx(okp) - xt(:,4), -P.vy(okp) - xt(:,5), -P.vz(okp) - xt(:,6)];
X.est_rate_hz = 1/median(diff(tp));
X.est_pos_rms = sqrt(mean(ep.^2,1)); X.est_vel_rms = sqrt(mean(ev.^2,1));
X.est_pos_sd = std(ep,1,1); X.est_vel_sd = std(ev,1,1);
tq = seconds(Q.timestamp) - L.tArm; okq = tq >= win(1) & tq < min(win(2), tL); tq = tq(okq); q = double(Q.q(okq,:));
eul = [atan2(2*(q(:,1).*q(:,2)+q(:,3).*q(:,4)), 1-2*(q(:,2).^2+q(:,3).^2)), ...
       -asin(max(-1,min(1,2*(q(:,1).*q(:,3)-q(:,4).*q(:,2))))), ...
       -atan2(2*(q(:,1).*q(:,4)+q(:,2).*q(:,3)), 1-2*(q(:,3).^2+q(:,4).^2))];
xa = interp1(L.t, L.x(:,7:9), tq); ea = eul - xa; ea(:,3) = angle(exp(1i*ea(:,3)));
X.est_att_rms = sqrt(mean(ea.^2,1)); X.est_att_sd = std(ea,1,1);
X.est_att_mean = mean(ea,1); X.est_pos_mean = mean(ep,1);
% effective delay of the estimate (min mean-square difference over lags 0-40 ms)
lags = 0:40; cp = zeros(size(lags)); ca = cp;
for k = 1:numel(lags)
    xs = interp1(L.t, L.x(:,1:2), tp - lags(k)*1e-3); cp(k) = mean(sum(([P.x(okp), -P.y(okp)] - xs).^2,2));
    as = interp1(L.t, L.x(:,7:8), tq - lags(k)*1e-3); ca(k) = mean(sum((eul(:,1:2) - as).^2,2));
end
[~,k1] = min(cp); [~,k2] = min(ca); X.lag_pos_ms = lags(k1); X.lag_att_ms = lags(k2);
% PSD of the x-position error and roll-angle error
fs = round(X.est_rate_hz); tg = (win(1):1/fs:win(2)-1/fs)';    % resampled at the logged estimate rate
exg = interp1(tp, ep(:,1), tg, 'linear', 'extrap'); erg = interp1(tq, ea(:,1), tg, 'linear', 'extrap');
[X.psd_pos, X.psd_f] = pwelch(exg - mean(exg), 512, [], [], fs);
X.psd_att = pwelch(erg - mean(erg), 512, [], [], fs);
end
