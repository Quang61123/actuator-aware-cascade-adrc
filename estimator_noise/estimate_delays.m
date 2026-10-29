function T = estimate_delays()
%ESTIMATE_DELAYS Effective delay of each estimated signal.
%   S6 runs with EKF2 input (est_s6_*, PID-Budget, AA-Budget, AA-HalfGain,
%   seeds 91001-91003), window 10-40 s: for every signal j the delay d
%   (0-20 ms) minimizing the mean-square difference between the signal used
%   by the controller and the delayed true signal, i.e. the lag of the best
%   delayed-truth approximation of the estimate. Hover runs are not used
%   because the true motion is too small to identify a lag.
%   Writes ../results/estimator_noise/estimate_delays.csv (per run) and
%   est_delays_summary.csv (median over runs).
here = fileparts(mfilename('fullpath'));
DW = fullfile(px4sitl_work_dir(),'est_data');
RES = fullfile(here,'..','results','estimator_noise');
f = dir(fullfile(DW,'est_s6_*.mat'));
names = {'x','y','z','vx','vy','vz','phi','theta','psi','phidot','thetadot','psidot'};
lags = 0:20; T = table();
for r = 1:numel(f)
    best = zeros(1,12); mse = zeros(numel(lags),12);
    for a = 1:numel(lags)
        E = estimation_errors(fullfile(DW,f(r).name), lags(a)*ones(1,12), [10 40]);
        mse(a,:) = mean(E.e.^2,1);
    end
    [~,ib] = min(mse,[],1); best(:) = lags(ib);
    rel = sqrt(min(mse,[],1)./mse(lags==0,:));                 % error reduction vs. zero lag
    T = [T; [table(string(f(r).name), E.meta.ctrl, E.meta.seed, 'VariableNames', {'Run','Ctrl','Seed'}), ...
        array2table(best, 'VariableNames', strcat('d_', names)), array2table(rel, 'VariableNames', strcat('rmsratio_', names))]]; %#ok<AGROW>
end
writetable(T, fullfile(RES,'estimate_delays.csv'));
Sm = array2table(median(T{:, strcat('d_', names)},1), 'VariableNames', strcat('d_', names));
writetable(Sm, fullfile(RES,'est_delays_summary.csv'));
end
