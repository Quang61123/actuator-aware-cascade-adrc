function csn = chan_update_base(kind, gj, cs, y, v, rr, ua, Ts, ilim)
%CHAN_UPDATE_BASE One-sample update of the channel state using the APPLIED
%(post-saturation) virtual input ua.  Forward-Euler discretisation.
%#codegen
csn = cs;
switch kind
    case {1,2}
        i = cs(1) + Ts*(rr - y);
        if gj(2) > 0
            % anti-windup clamp: steady-state share required by the set-point
            % weight, kp*(1-bsp)*|r|, plus a disturbance allowance ilim
            lim = (ilim + gj(1)*(1 - gj(5))*abs(rr))/gj(2);
            i = min(max(i, -lim), lim);
        end
        csn(1) = i;
        if gj(4) > 0
            a = exp(-gj(4)*Ts); csn(2) = a*cs(2) + (1-a)*v;
        end
    case 3
        wo = gj(2); b1 = 3*wo; b2 = 3*wo^2; b3 = wo^3;
        e = y - cs(1);
        csn(1) = cs(1) + Ts*(cs(2) + b1*e);
        csn(2) = cs(2) + Ts*(cs(3) + ua + b2*e);
        csn(3) = cs(3) + Ts*(b3*e);
    otherwise
        wo = gj(2); b1 = 2*wo; b2 = wo^2;
        e = v - cs(1);
        csn(1) = cs(1) + Ts*(cs(2) + ua + b1*e);
        csn(2) = cs(2) + Ts*(b2*e);
end
end


