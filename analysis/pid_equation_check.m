function T = pid_equation_check()
%PID_EQUATION_CHECK Sample-by-sample check of the PID baseline equations of
%   the manuscript (Section II, "PID baseline", Eqs. (pid)-(pidstate))
%   against the production implementation (chan_output.m, kind 2, and
%   chan_update_base.m, kind 2) that ctrl_law calls in MATLAB, Simulink and
%   the code generated for PX4.
%
%   Manuscript equations (implemented independently below):
%     u_k        = Kp*(b*r_k - y_k) + Ki*iota_k + Kd*(rd_k - vf_k) + ff*rdd_k
%     iota_{k+1} = sat_{l_k}( iota_k + Ts*(r_k - y_k) ),
%                  l_k = (ilim + Kp*(1-b)*|r_k|)/Ki
%     vf_{k+1}   = af*vf_k + (1-af)*v_k,  af = exp(-wf*Ts)
%     iota_0 = 0, vf_0 = 0  (all controller states start at zero, as in
%     run_quad_sim/run_quad_hover/run_quad_stress, the Simulink model and
%     sitl_ctrl_step, where they are also reset while disarmed).
%   ilim = 2, 3, 20, 5 for (xy, z, phi-theta, psi) as in ctrl_law.
%   Attitude classes use rd = rdd = 0 and ff = 0, as in ctrl_law.
%
%   Gains: results/analysis/setpoint_designs.mat, C(k).star (PID-Budget (b*),
%   PID-Budget-S, PID-E; full precision).
%   Output: results/analysis/pid_equation_consistency.csv
%
%   The test sequence (seed 1, 20000 steps) combines a sinusoid, steps and a
%   large transient error so that the anti-windup clamp is active.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'MATLAB_simulation'));
A = fullfile(root, 'results', 'analysis');
S = load(fullfile(A, 'setpoint_designs.mat'));
cls = {'x', 1, 2; 'z', 3, 3; 'phi', 4, 20; 'psi', 6, 5};   % class, channel index, ilim
Ts = 1e-3; N = 20000; rng(1);
rows = {};
for d = 1:numel(S.C)
    G = S.C(d).star.G; nm = S.C(d).star.name;
    assert(S.C(d).star.kind == 2, 'PID design expected (kind 2)');
    for c = 1:size(cls, 1)
        j = cls{c, 2}; gj = G(6*(j-1) + (1:6)); il = cls{c, 3}; att = j >= 4;
        kp = gj(1); ki = gj(2); kd = gj(3); wf = gj(4); bb = gj(5); ff = G(37)*(~att);
        assert(wf > 0 && ki > 0, 'filtered-rate PID with integral action expected');
        t = (0:N-1)'*Ts;
        r = 2*sin(0.7*t) + 1.5*(t > 5) - 3*(t > 12);
        rd = 1.4*cos(0.7*t); rdd = -0.98*sin(0.7*t);
        if att, rd = 0*rd; rdd = 0*rdd; end
        y = r + 0.5*randn(N,1).*(t < 3) + 0.05*randn(N,1) + 0.8*(t > 8 & t < 9);
        v = randn(N,1);
        cs = zeros(4,1);            % production state: cs(1) = iota, cs(2) = vf
        io = 0; vf = 0;             % manuscript state: iota_0 = 0, vf_0 = 0
        eu = 0; eI = 0; ev = 0; ncl = 0;
        for k = 1:N
            u  = chan_output(2, gj, cs, y(k), v(k), r(k), rd(k), rdd(k), ff);
            up = kp*(bb*r(k) - y(k)) + ki*io + kd*(rd(k) - vf) + ff*rdd(k);
            eu = max(eu, abs(u - up));
            csn = chan_update_base(2, gj, cs, y(k), v(k), r(k), 0, Ts, il);
            lim = (il + kp*(1 - bb)*abs(r(k)))/ki;
            iraw = io + Ts*(r(k) - y(k)); ion = min(max(iraw, -lim), lim);
            ncl = ncl + (ion ~= iraw);
            af = exp(-wf*Ts); vfn = af*vf + (1 - af)*v(k);
            eI = max(eI, abs(csn(1) - ion)); ev = max(ev, abs(csn(2) - vfn));
            cs = csn; io = ion; vf = vfn;
        end
        rows(end+1, :) = {nm, cls{c,1}, kp, ki, kd, wf, bb, ff, il, N, ncl, eu, eI, ev}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'design','class','Kp','Ki','Kd','wf','b','ff','ilim', ...
    'steps','clamp_active_steps','max_abs_err_u','max_abs_err_integrator','max_abs_err_ratefilter'});
writetable(T, fullfile(A, 'pid_equation_consistency.csv'));
fprintf('pid_equation_check: max |u| error %.3g, max |iota| error %.3g, max |vf| error %.3g (%d cases)\n', ...
    max(T.max_abs_err_u), max(T.max_abs_err_integrator), max(T.max_abs_err_ratefilter), height(T));
end
