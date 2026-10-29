function [g, ok, q] = pid_repair(g, c, S, gin, dubc, req, oa)
%PID_REPAIR Bandwidth and set-point-weight repair of a PID class gain
%   vector (general step of the PID re-searches; same procedure as
%   feasibility_repair, extended to the response and stack checks):
%     for j = 1..70: scale the bandwidth of the original gains by factor
%     (scale_design: Kp t^2, Ki t^3, Kd t, wf t); if OS > 10 %, reduce the
%     set-point weight b by 0.8 up to 8 times (b changes the reference
%     response only) and keep it; accept when sigma_u <= ub(c),
%     sigma_du <= dubc, Ms <= 1.98 (4000 frequencies plus refinement),
%     OS <= 10 % and 2% settling <= req(c) (400 s), Schur stable;
%     factor *= 0.985 for the first 10 steps, then *= 0.9.
%   ok = false if no step passes. q = [sigma_e sigma_u sigma_du Ms OS ts rho].
orig = g; factor = 1; ok = false; q = [];
for j = 1:70
    g = scale_design(2, orig, factor);
    [os, ts] = step_response_metrics(chan_model(2, g, S.Ts, gin, S.act, [], oa), S.Ts);
    if isfinite(os) && os > S.osmax
        for ib = 1:8
            g(5) = g(5)*0.8; [os, ts] = step_response_metrics(chan_model(2, g, S.Ts, gin, S.act, [], oa), S.Ts);
            if os <= S.osmax, break, end
        end
        orig(5) = g(5);
    end
    z = channel_stats(2, g, c, S, gin, oa, S.act, false, 4000);
    q = [z(1:4) os ts z(6)];
    if all(isfinite(q)) && z(2) <= S.ub(c)*(1+1e-9) && z(3) <= dubc*(1+1e-9) && z(4) <= S.Msmax*(1+1e-6) ...
            && os <= S.osmax && ts <= req(c) && z(6) < 1
        ok = true; return
    end
    if j <= 10, factor = factor*0.985; else, factor = factor*0.9; end
end
end
