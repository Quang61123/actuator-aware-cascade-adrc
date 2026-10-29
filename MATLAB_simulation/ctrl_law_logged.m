function [T, tau, att_ref, csn, ulog] = ctrl_law_logged(kind, G, cs, meas, r, pv, Ts)
%CTRL_LAW_LOGGED Cascade controller (structure and inversion).
%   kind 1 PID (reference gains), 2 optimized PID,
%        4 RO-LADRC (reduced-order LADRC, ESO input = commanded control),
%        7 AA RO-LADRC (actuator-aware, delay-aligned inner/z observers;
%          attitude-aided outer observer).
%        8 HG-RO half-gain predictive observer; 9 HG-AA half-gain actuator-aware.
%   G : 37x1 (6 gains per channel, G(37) feed-forward flag)
%   cs: 96x1 (16 per channel).  meas: delayed, noisy measurement (12x1).
%#codegen
m = pv(1); g = pv(2); I = pv(7:9);
alpha = exp(-Ts/pv(16)); nd = round(pv(17));
ff = G(37); ilim = [2 2 3 20 20 5]; angmax = 0.35;
csn = cs; ua = zeros(6,1);
lk = kind; if kind == 7 || kind == 8 || kind == 9, lk = 4; end            % same output law RO-LADRC
ux = chan_output(lk, G(1:6),  cs(1:4),   meas(1), meas(4), r(1), r(5), r(8), ff);
uy = chan_output(lk, G(7:12), cs(17:20), meas(2), meas(5), r(2), r(6), r(9), ff);
psi = meas(9); cp = cos(psi); sp = sin(psi);
thd = min(max((ux*cp + uy*sp)/g, -angmax), angmax);
phd = min(max((ux*sp - uy*cp)/g, -angmax), angmax);
if kind == 4 || kind == 7 || kind == 8 || kind == 9      % attitude-aided outer ESO 
    ua(1) = g*(meas(8)*cp + meas(7)*sp);  ua(2) = g*(meas(8)*sp - meas(7)*cp);
else
    ua(1) = g*(thd*cp + phd*sp);  ua(2) = g*(thd*sp - phd*cp);
end
uz = chan_output(lk, G(13:18), cs(33:36), meas(3), meas(6), r(3), r(7), r(10), ff);
cc = cos(meas(7))*cos(meas(8));
T = min(max(m*(g + uz)/cc, 0), pv(12));
ua(3) = T*cc/m - g;
aref = [phd; thd; r(4)];
uat = zeros(3,1);
for j = 1:3
    ic = 48 + 16*(j-1) + (1:4); ig = 18 + 6*(j-1) + (1:6);
    uat(j) = chan_output(lk, G(ig), cs(ic), meas(6+j), meas(9+j), aref(j), 0, 0, 0);
end
tau = min(max(I(:).*uat, -pv(13:15)), pv(13:15));
ua(4:6) = tau./I(:);
att_ref = aref;
yv = [meas(1:3); meas(7:9)]; vv = [meas(4:6); meas(10:12)]; rv = [r(1:3); aref];
for j = 1:6
    ic = 16*(j-1) + (1:16); ig = 6*(j-1) + (1:6);
    kj = kind; if kind == 7 && j <= 2, kj = 4; end
    if kind==9 && j<=2,kj=8;end   % x-y: attitude-aided (no motor model)
    csn(ic) = chan_update(kj, G(ig), cs(ic), yv(j), vv(j), rv(j), ua(j), Ts, ilim(j), alpha, nd);
end
ulog=[g*(thd*cp+phd*sp);g*(thd*sp-phd*cp);ua(3:6)];
end


