function T = rotor_increment_crosscov()
%ROTOR_INCREMENT_CROSSCOV Effect of the cross-covariance between the virtual-control
%   increments on the rotor-command increment J.
%   J^2 = trace(Cov(dw_c))/4, dw_c = Jm * diag([m Ixx Iyy Izz]) * du_v,
%   Jm = linearized inverse mixer about hover (rad/s per N or N m).
%   From the logged virtual-control increments (debug_array, 1 kHz, hover
%   noise-only, window 20-60 s) of every hover run, J is computed (i) with
%   the full measured covariance of du_v = [du_z du_phi du_theta du_psi] and
%   (ii) with its diagonal only (the uncorrelated-class assumption of the
%   prediction); both are compared with the logged unquantized J.
%   Data: est_data/est_hover_* and redesign_data/E0H_*.
%   Writes ../results/analysis/rotor_increment_crosscov.csv.
here = fileparts(mfilename('fullpath')); RA = fullfile(here,'..','results','analysis');
pv = quad_params(); k = pv(4); l = pv(3); b = pv(5); wh = sqrt(pv(1)*pv(2)/(4*k));
Mx = [k k k k; 0 -l*k 0 l*k; -l*k 0 l*k 0; -b b -b b]; Jm = diag(1/(2*wh)*ones(4,1))/Mx;
Sc = diag([pv(1) pv(7:9)']);
f = [dir(fullfile(px4sitl_work_dir(),'est_data','est_hover_*.mat')); dir(fullfile(px4sitl_work_dir(),'redesign_data','E0H_hover_*.mat'))];
T = table();
for r = 1:numel(f)
    D = load(fullfile(f(r).folder, f(r).name), 'R', 'c', 'seed'); R = D.R;
    kw = R.tc >= 20 & R.tc < 60; okw = kw & R.ok; okw(1) = false; ip = find(okw);
    dU = double(R.ua(ip,:) - R.ua(ip-1,:)); dU = dU(:, [4 1 2 3]);          % [z phi theta psi]
    C = cov(dU, 1); A = Jm*Sc;
    Jfull = sqrt(trace(A*C*A')/4); Jdiag = sqrt(trace(A*diag(diag(C))*A')/4);
    dW = pv(11)*double(R.wn(ip,:) - R.wn(ip-1,:)); Jlog = sqrt(mean(dW(:).^2));
    T = [T; table(D.c, D.seed, Jlog, Jfull, Jdiag, 'VariableNames', {'Ctrl','Seed','J_logged','J_fullCov','J_diagCov'})]; %#ok<AGROW>
end
T.DiagOverFull = T.J_diagCov ./ T.J_fullCov; T.FullOverLogged = T.J_fullCov ./ T.J_logged;
writetable(T, fullfile(RA,'rotor_increment_crosscov.csv'));
end
