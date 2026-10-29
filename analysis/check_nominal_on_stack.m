function V = check_nominal_on_stack()
%CHECK_NOMINAL_ON_STACK Checks the nominal designs PID-Budget, AA-Budget and AA-HalfGain
%   against the estimator-stack specification used for the -E designs (verify_design_stack: stack delay convention,
%   Ms <= 2 / target 1.98, OS, settling, local white ceilings, frozen Bcas and Bpsi with margin 1),
%   with the same pooled noise model as the -E design. Writes check_nominal_on_stack.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub;
L = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF = L.NF;
F = load(fullfile(RA,'frozen_ceilings_and_delays.mat'));
V = table();
for i = [1 2 4], co = D.C(i); V = [V; verify_design_stack(co, S, dub, NF, F.ds, F.Bcas, F.Bpsi, 4000, 1)]; end %#ok<AGROW>
writetable(V, fullfile(RA,'check_nominal_on_stack.csv')); disp(V(:,{'Design','Class','Ms_E','OS_E','ts_E','sigma_du_est','Ceiling_est','Pass','MsTarget198'}))
end