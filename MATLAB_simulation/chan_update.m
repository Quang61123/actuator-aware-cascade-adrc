function csn = chan_update(kind, gj, cs, y, v, rr, ua, Ts, ilim, alpha, nd)
%CHAN_UPDATE One-sample update of a channel state (16x1).
%   kind 1,2 : PID        kind 4 : RO-LADRC (ESO input = ua)
%   kind 7   : actuator-aware, delay-aligned RO-LADRC
%              cs(1:2) ESO [rate est., disturbance est.]
%              cs(3)   motor-model output a_hat(k) (first-order, pole alpha)
%              cs(4:13) history a_hat(k-1..k-10) for delay alignment
%   ESO input for kind 7 is a_hat(k-nd): the acceleration actually applied
%   at the instant the (delayed) measurement refers to.
%#codegen
csn = cs;
switch kind
    case 8
        L=observer_half_gain(gj(2),Ts);e=v-cs(1);
        csn(1)=cs(1)+Ts*(cs(2)+ua)+L(1)*e;
        csn(2)=cs(2)+L(2)*e;
    case {1,2,3,4}
        csn(1:4) = chan_update_base(kind, gj, cs(1:4), y, v, rr, ua, Ts, ilim);
    otherwise
        wo = gj(2); b1 = 2*wo; b2 = wo^2;
        if kind==9,L=observer_half_gain(wo,Ts);b1=L(1)/Ts;b2=L(2)/Ts;end
        m = cs(3);
        if nd >= 1, uin = cs(3 + nd); else, uin = m; end
        e = v - cs(1);
        csn(1) = cs(1) + Ts*(cs(2) + uin + b1*e);
        csn(2) = cs(2) + Ts*(b2*e);
        csn(5:13) = cs(4:12);
        csn(4) = m;
        csn(3) = alpha*m + (1 - alpha)*ua;
end
end


