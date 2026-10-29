function [DE, V] = finalize_estimator_designs()
%FINALIZE_ESTIMATOR_DESIGNS Assemble and finalize the estimator-stack (-E) designs.
%   (i) The AA-Budget-E and AA-HalfGain-E designs are taken from the run
%       saved in design_E_checkpoint_AA_backup.mat (the resumed run that
%       redid PID-E reloaded the unmodified base designs for these two
%       entries; the searches themselves are unchanged).
%   (ii) PID-E: the SQP solution meets OS <= 10 % only to the solver
%       tolerance (10.000001 %). As in the repair step of the primary design
%       (feasibility_repair), only the set-point weight b of the affected PID
%       classes is reduced in steps of 0.1 % until OS <= 10 %; b changes the
%       reference response only, not the noise covariances or Ms.
%       For the PID-E yaw class neither b nor bandwidth scaling lowered OS
%       below 10 %; its search was repeated with a 9.99 % OS design target
%       (pid_estimator_yaw_margin) and that result is used.
%   Re-verifies all designs (verify_design_stack, 4000 frequencies) and writes
%   designs_E.mat (DE, V, ...) and designs_E_verify.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis');
D = load(fullfile(RP,'design_gains.mat')); S = D.S; dub = D.dub;
L = load(fullfile(RP,'estimator_noise','est_noise_models.mat')); NF = L.NF;
% raw search outputs: PID-E from the checkpoint of its run, AA designs from their run
P = load(fullfile(RA,'design_E_checkpoint.mat')); A = load(fullfile(RA,'design_E_checkpoint_AA_backup.mat'));
F = load(fullfile(RA,'frozen_ceilings_and_delays.mat'));
DE = P.DE; DE(2) = A.DE(2); DE(3) = A.DE(3); ds = F.ds; Bcas = F.Bcas; Bpsi = F.Bpsi;
K.search = P.search(startsWith(P.search.Design, "PID-E"), :);
fy = fullfile(RA,'pid_estimator_yaw_margin.mat');   % yaw search repeated with a 9.99 % OS target (pid_estimator_yaw_margin)
if exist(fy,'file'), Y = load(fy); if ~isempty(Y.best), DE(1).G(31:36) = Y.best'; end, end
co = DE(1); oa = [co.tauhat co.ndhat]; log = table();
for h = [1 3 4]
    ix = {1:6, [], 19:24, 31:36}; ii = ix{h}; gin = []; if h == 1, gin = co.G(19:24)'; end
    g = co.G(ii)'; b0 = g(5);
    r = est_stack_tools('resp', 2, g, h, gin, oa, ds); n = 0;
    while r(1) > S.osmax && n < 200, g(5) = g(5)*0.999; n = n + 1; r = est_stack_tools('resp', 2, g, h, gin, oa, ds); end
    co.G(ii) = g'; if h == 1, co.G(7:12) = g'; elseif h == 3, co.G(25:30) = g'; end
    log = [log; table(h, b0, g(5), n, r(1), r(2), 'VariableNames', {'Class','b_before','b_after','Steps','OS_E','ts_E'})]; %#ok<AGROW>
end
DE(1) = co;
V = verify_design_stack(D.C(1), S, dub, NF, ds, Bcas, Bpsi, 4000, 1); V.Design(:) = "PID-Budget (reference)";
for j = 1:3, V = [V; verify_design_stack(DE(j), S, dub, NF, ds, Bcas, Bpsi, 4000)]; end %#ok<AGROW>
search = [K.search; A.search];
save(fullfile(RA,'designs_E.mat'), 'DE', 'V', 'search', 'Bcas', 'Bpsi', 'ds', 'log');
writetable(V, fullfile(RA,'designs_E_verify.csv')); writetable(search, fullfile(RA,'design_E_search.csv'));
writetable(log, fullfile(RA,'pidE_setpoint_weight.csv'));
end
