function M = chan_model(kind, gj, Ts, gin, act, bs, oact)
%CHAN_MODEL Implemented sampled-data linear model of one normalised channel
%   including first-order actuator lag and a measurement delay, built from
%   the SAME controller code (chan_output / chan_update).
%   act = [tau_m, nd] of the PLANT; bs = true/nominal input-gain ratio (default 1);
%   oact = [tau_m, nd] assumed INSIDE the controller (actuator-aware ESO of
%   kind 7); default oact = act. Use oact = design values for off-nominal
%   robustness analysis (the observer is not retuned).
%   Single loop : y'' = d + a,  a(k+1) = alpha a(k) + (1-alpha) bs u_t(k)
%   Cascade     : y'' = d + g*theta, inner theta loop with its own lag,
%                 delay and controller (gains gin); outer ESO of kinds 4/7
%                 is attitude-aided (input g*theta_measured).
%   Inputs  w = [w_d, n_y, n_v, n_th, n_om, r, w_u]; outputs y, u, u_t.
if nargin < 4, gin = []; end
if nargin < 6 || isempty(bs), bs = 1; end
if nargin < 7 || isempty(oact), oact = act; end
g = 9.81; casc = ~isempty(gin);
alpha = exp(-Ts/act(1)); nd = round(act(2));
alo = exp(-Ts/oact(1)); ndo = round(oact(2));
lk = kind; if kind == 7 || kind == 8 || kind == 9 || kind == 10, lk = 4; end
ko = kind; if casc && kind == 7, ko = 4; end
if casc && kind==9,ko=8;end          % outer x-y law
cmdin = false; if casc && kind == 10, ko = 4; cmdin = true; end   % A4-CommandedTilt, commanded x-y ESO input
nx1 = 4 + 2*nd + 16; nx = nx1 + casc*(3 + 2*nd + 16); nw = 7;
A = zeros(nx); B = zeros(nx,nw); Cy = zeros(1,nx); Cu = zeros(1,nx); Ct = zeros(1,nx);
Du = zeros(1,nw); Dt = zeros(1,nw);
for c = 1:nx, e = zeros(nx,1); e(c) = 1; [A(:,c), Cy(c), Cu(c), Ct(c)] = f1(e, zeros(nw,1)); end
for c = 1:nw, e = zeros(nw,1); e(c) = 1; [B(:,c), ~, Du(c), Dt(c)] = f1(zeros(nx,1), e); end
% remove decoupled, unused states (identity rows/cols without input)
keep = true(1,nx); keep(3) = false;          % d (random walk) is handled separately
for i = 4:nx
    ei = zeros(1,nx); ei(i) = 1;
    if isequal(A(i,:), ei) && isequal(A(:,i)', ei) && ~any(B(i,:)), keep(i) = false; end
end
M.A = A; M.B = B; M.Cy = Cy; M.Cu = Cu; M.Ct = Ct; M.Du = Du; M.Dt = Dt;
M.keep = find(keep); M.id = 3;
    function [Xn, yo, uo, to] = f1(X, w)
        y = X(1); v = X(2); d = X(3); a = X(4);
        yb = X(5:4+nd); vb = X(5+nd:4+2*nd); cs = X(5+2*nd:nx1);
        if nd >= 1, ym = yb(nd) + w(2); vm = vb(nd) + w(3); else, ym = y + w(2); vm = v + w(3); end
        u  = chan_output(lk, gj, cs(1:4), ym, vm, w(6), 0, 0, 0);
        ut = u + w(7);
        if ~casc
            csn = chan_update(kind, gj, cs, ym, vm, w(6), u, Ts, inf, alo, ndo);
            yn = y + Ts*v + Ts^2/2*(d + a);  vn = v + Ts*(d + a);
            an = alpha*a + (1-alpha)*bs*ut;  Xi = [];
        else
            o = nx1;
            th = X(o+1); om = X(o+2); ath = X(o+3);
            thb = X(o+4:o+3+nd); omb = X(o+4+nd:o+3+2*nd); ci = X(o+4+2*nd:nx);
            if nd >= 1, thm = thb(nd) + w(4); omm = omb(nd) + w(5); else, thm = th + w(4); omm = om + w(5); end
            if (ko == 4 || ko == 8) && ~cmdin, uae = g*thm; else, uae = u; end
            csn = chan_update(ko, gj, cs, ym, vm, w(6), uae, Ts, inf, alo, ndo);
            thd = ut/g;
            ui  = chan_output(lk, gin, ci(1:4), thm, omm, thd, 0, 0, 0);
            cin = chan_update(kind, gin, ci, thm, omm, thd, ui, Ts, inf, alo, ndo);
            yn = y + Ts*v + Ts^2/2*(d + g*th);  vn = v + Ts*(d + g*th);  an = a;
            thn = th + Ts*om + Ts^2/2*ath;  omn = om + Ts*ath;  athn = alpha*ath + (1-alpha)*bs*ui;
            Xi = [thn; omn; athn; shift(thb, th); shift(omb, om); cin];
        end
        Xn = [yn; vn; d + w(1); an; shift(yb, y); shift(vb, v); csn; Xi];
        yo = y; uo = u; to = ut;
    end
end
function b = shift(b, newest)
if isempty(b), return, end
b = [newest; b(1:end-1)];
end



