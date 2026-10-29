function u = chan_output(kind, gj, cs, y, v, rr, rd, rdd, ff)
%CHAN_OUTPUT Virtual-acceleration command of one channel  y'' = d + u.
%   kind 1,2 : 2-DOF PID (gj = [kp ki kd wf bsp]); D-term on (filtered) rate
%              measurement, set-point weight bsp on the proportional term
%   kind 3   : full-order LADRC (gj = [wc wo]), ESO driven by position only
%   kind 4   : rate-aided reduced-order LADRC (gj = [wc wo])
%   cs = 4x1 channel state, y/v = measured position/rate, rr/rd/rdd = reference.
%#codegen
switch kind
    case {1,2}
        if gj(4) > 0, vf = cs(2); else, vf = v; end
        u = gj(1)*(gj(5)*rr - y) + gj(2)*cs(1) + gj(3)*(rd - vf) + ff*rdd;
    case 3
        kp = gj(1)^2; kd = 2*gj(1);
        u = kp*(rr - cs(1)) + kd*(rd - cs(2)) + ff*rdd - cs(3);
    otherwise
        kp = gj(1)^2; kd = 2*gj(1);
        u = kp*(rr - y) + kd*(rd - cs(1)) + ff*rdd - cs(2);
end
end


