function V = verify_design(co, S, dub, req, cap)
%VERIFY_DESIGN Full check of a design against the covariance design constraints plus the response
%   requirement, on the white-noise design model (design delays):
%   per class sigma_e, sigma_u, local sigma_du, Ms (4000 frequencies plus
%   refinement), spectral radius, OS and 2% settling (400 s; xy on the
%   cascade), and the cascade inner increment against cap.
oa = [co.tauhat co.ndhat]; nm = {'xy','z','phth','psi'};
gs = {co.G(1:6)', co.G(13:18)', co.G(19:24)', co.G(31:36)'}; gin = {co.G(19:24)', [], [], []};
zc = cascade_stats(co, S); V = table();
for h = 1:4
    q = channel_stats(co.kind, gs{h}, h, S, gin{h}, oa, S.act, false, 4000);
    [os, ts] = step_response_metrics(chan_model(co.kind, gs{h}, S.Ts, gin{h}, S.act, [], oa), S.Ts);
    cdu = NaN; okc = true; if h == 1, cdu = zc(4); okc = cdu <= cap*(1-0.01) + 1e-12; end
    % final check: Ms <= 2 (design target 1.98, reported in MsTarget198 to the solver tolerance 1e-6)
    pass = all(isfinite(q(1:4))) && q(2) <= S.ub(h)*(1+1e-9) && q(3) <= dub(h)*(1+1e-9) && q(4) <= 2 && q(6) < 1 ...
        && os <= S.osmax && ts <= req(h) && okc;
    V = [V; table(string(co.name), string(nm{h}), gs{h}(1), gs{h}(2), q(1), q(2), q(3), q(4), q(6), os, ts, cdu, pass, q(4) <= S.Msmax*(1+1e-6), ...
        'VariableNames', {'Design','Class','g1','g2','sigma_e','sigma_u','sigma_du','Ms','rho','OS','ts','cascade_du','Pass','MsTarget198'})]; %#ok<AGROW>
end
end
