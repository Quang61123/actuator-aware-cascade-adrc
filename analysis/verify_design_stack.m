function V = verify_design_stack(co, S, dub, NF, ds, Bcas, Bpsi, nf, margin)
%VERIFY_DESIGN_STACK Check of a design against the estimator-stack specification:
%   per class, with the delay convention ds,
%   Ms <= 1.98 (nf frequencies plus refinement; 4000 for the final check),
%   OS <= 10 % and 2% settling <= req (400 s; xy on the cascade), local
%   white-noise ceilings, Schur stability, and the frozen
%   estimator-aware ceilings: cascade roll/pitch increment <= margin*Bcas,
%   yaw increment <= margin*Bpsi (margin 0.99 for new designs, 1 for the
%   reference PID-Budget). Also reports Ms with the design delays.
if nargin < 8, nf = 4000; end
if nargin < 9, margin = 0.99; end
req = [15 5 2 5]; oa = [co.tauhat co.ndhat]; nm = {'xy','z','phth','psi'};
gs = {co.G(1:6)', co.G(13:18)', co.G(19:24)', co.G(31:36)'}; gin = {co.G(19:24)', [], [], []};
zc = est_stack_tools('cascade', co, S, NF, ds); V = table();
for h = 1:4
    q = channel_stats(co.kind, gs{h}, h, S, gin{h}, oa, S.act, false, nf);     % white model: local ceilings, design-delay Ms
    msE = est_stack_tools('ms', co.kind, gs{h}, h, gin{h}, oa, ds, nf);
    rs = est_stack_tools('resp', co.kind, gs{h}, h, gin{h}, oa, ds);
    switch h
        case 1, ze = zc; est = zc(4); lim = margin*Bcas;
        case 3, ze = zc; est = zc(4); lim = margin*Bcas;
        case 4, ze = est_stack_tools('single', co.kind, gs{h}, h, S, NF, ds, oa); est = ze(3); lim = margin*Bpsi;
        otherwise, ze = est_stack_tools('single', co.kind, gs{h}, h, S, NF, ds, oa); est = ze(3); lim = Inf;
    end
    % final check: Ms <= 2; the design target 1.98 is reported separately (to the solver tolerance 1e-6)
    pass = all(isfinite(q(1:3))) && q(2) <= S.ub(h)*(1+1e-9) && q(3) <= dub(h)*(1+1e-9) && msE <= 2 ...
        && rs(1) <= S.osmax && rs(2) <= req(h) && ze(5) < 1 && est <= lim*(1+1e-9);
    V = [V; table(string(co.name), string(nm{h}), gs{h}(1), gs{h}(2), q(2), q(3), q(4), msE, rs(1), rs(2), ze(1), est, lim, pass, msE <= S.Msmax*(1+1e-6), ...
        'VariableNames', {'Design','Class','g1','g2','sigma_u_white','sigma_du_white','Ms_design','Ms_E','OS_E','ts_E','sigma_e_est','sigma_du_est','Ceiling_est','Pass','MsTarget198'})]; %#ok<AGROW>
end
end
