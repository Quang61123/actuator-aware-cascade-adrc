function M = sitl_metrics(L)
%SITL_METRICS Metrics of one SITL run with the MATLAB definitions
%   (run_metrics: S6 window 10-40 s, S5 window 5-45 s). Virtual-input
%   clipping is internal to PX4 and is not available (NaN).
pv = quad_params();
if L.scen == 5, win = [5 45]; else, win = [10 40]; end
M = run_metrics(L, win, pv);
M.win = win;
end
