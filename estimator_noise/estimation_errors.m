function E = estimation_errors(file, dsig, win)
%ESTIMATION_ERRORS Estimator errors at controller rate for one SITL run with controller-rate logging.
%   file : <work dir>\est_data\<tag>.mat (run_sitl_estimator)
%   dsig : 1x12 delays [samples] of the 12 measured signals (x y z vx vy vz
%          phi theta psi phidot thetadot psidot); the error of signal j is
%          e_j(k) = meas_j(k) - truth_j(k - dsig(j)),
%          i.e. the measurement is modelled as a delayed true signal plus
%          noise, as in the covariance model (white model: dsig = 5).
%   win  : window [s] of controller time.
%   Controller step k (time k*Ts) and plant sample k+1 refer to the same
%   instant (verified with the estimator-bypass calibration run, in which
%   meas(k) = x(k+1-5) + noise exactly).
%   Truth: model-frame position, velocity, Euler angles and Euler-angle
%   rates W^-1(eta) * body rates of the plant state.
%   E.e (N x 12) errors, E.tc times, E.seg segment index (consecutive steps),
%   E.meta (ctrl, seed, hover, raw).
D = load(file);
R = D.R; x = double(D.Lx);
k = round(R.tc*1000); i = k + 1;
sel = R.tc >= win(1) & R.tc < win(2);
k = k(sel); i = i(sel); meas = double(R.meas(sel,:));
% Euler-angle rates of the plant
phi = x(:,7); th = x(:,8); p = x(:,10); q = x(:,11); r = x(:,12);
ed = [p + sin(phi).*tan(th).*q + cos(phi).*tan(th).*r, cos(phi).*q - sin(phi).*r, (sin(phi).*q + cos(phi).*r)./cos(th)];
truth = [x(:,1:9), ed];
e = zeros(numel(i),12);
for j = 1:12
    ii = i - dsig(j);
    e(:,j) = meas(:,j) - truth(ii,j);
end
E.e = e; E.tc = k/1000; E.seg = cumsum([1; diff(k) ~= 1]);
E.meta = struct('ctrl', D.c, 'seed', D.seed, 'hover', D.hv, 'raw', D.raw);
end
