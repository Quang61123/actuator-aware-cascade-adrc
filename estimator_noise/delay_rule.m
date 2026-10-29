function dsig = delay_rule()
%DELAY_RULE Signal delays used in the estimator-aware model.
%   From estimate_delays.csv (S6 runs): for each of the 12 signals the median over
%   the runs of the best lag and of the RMS ratio (error RMS at the best lag
%   / error RMS at zero lag). Rule: a signal is modelled as delayed by the
%   median lag only if this delay explains a material part of the error,
%   i.e. the median RMS ratio is at most 0.9; otherwise the estimate is
%   modelled without delay (meas = truth + coloured noise).
%   Writes est_delays_summary.csv (median lag, median ratio, chosen delay).
here = fileparts(mfilename('fullpath'));
RES = fullfile(here,'..','results','estimator_noise');
T = readtable(fullfile(RES,'estimate_delays.csv'));
names = {'x','y','z','vx','vy','vz','phi','theta','psi','phidot','thetadot','psidot'};
dmed = median(T{:, strcat('d_', names)}, 1);
rmed = median(T{:, strcat('rmsratio_', names)}, 1);
dsig = dmed .* (rmed <= 0.9);
Sm = table(string(names(:)), dmed(:), rmed(:), dsig(:), 'VariableNames', {'Signal','MedianLag_ms','MedianRMSRatio','ChosenDelay_ms'});
writetable(Sm, fullfile(RES,'est_delays_summary.csv'));
end
