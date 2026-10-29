function T = settling_check()
%SETTLING_CHECK Response check of all existing designs.
%   For every design (primary C(1:7), strict, -E, -E2), the unit reference
%   step of the implemented sampled linear channel model (chan_step, design
%   delays, white-noise model) of each class xy, z, phi/theta, psi gives
%   the 2% settling time and the overshoot. The horizon is 400 s, so that
%   slow loops are not cut off (the design search used 4*Tstep, which can
%   miss a late overshoot). Requirements: ts_xy <= 15 s, ts_z <= 5 s, ts_phth <= 2 s,
%   ts_psi <= 5 s; OS <= 10 %.
%   Writes ../results/analysis/settling_check.csv.
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
RP = fullfile(here,'..','results'); RA = fullfile(RP,'analysis'); if ~exist(RA,'dir'), mkdir(RA); end
D = load(fullfile(RP,'design_gains.mat')); S = D.S;
E = load(fullfile(RP,'estimator_noise','est_designs_E.mat')); E2 = load(fullfile(RP,'yaw_ceiling','yaw_ceiling_designs.mat'));
des = [D.C(1:7), D.Strict(1:3), E.DE(1:2), E2.DE2(1:2)];
for j = 1:numel(D.Curve)       % budget-curve candidates
    c = D.Curve(j).Controller; c.name = sprintf('Curve %s gamma=%.2f', D.Curve(j).Family, D.Curve(j).Gamma);
    des(end+1) = orderfields(rmfield(c, setdiff(fieldnames(c), fieldnames(des))), des(1)); %#ok<AGROW>
end
req = [15 5 2 5]; osmax = S.osmax; Tend = 400; T = table(); cls = {'xy','z','phth','psi'};
for r = 1:numel(des)
    co = des(r); oa = [co.tauhat co.ndhat];
    gs = {co.G(1:6)', co.G(13:18)', co.G(19:24)', co.G(31:36)'}; gin = {co.G(19:24)', [], [], []};
    ts = zeros(1,4); os = ts;
    for h = 1:4, [os(h), ts(h)] = step_response_metrics(chan_model(co.kind, gs{h}, S.Ts, gin{h}, S.act, [], oa), S.Ts, 1, Tend); end
    ok = ts <= req & os <= osmax + 1e-9;
    fails = string(strjoin(cls(~ok), ' '));
    T = [T; table(string(co.name), ts(1), ts(2), ts(3), ts(4), os(1), os(2), os(3), os(4), all(ok), fails, ...
        'VariableNames', {'Design','ts_xy','ts_z','ts_phth','ts_psi','OS_xy','OS_z','OS_phth','OS_psi','Pass','FailedClasses'})]; %#ok<AGROW>
end
writetable(T, fullfile(RA,'settling_check.csv'));
end
