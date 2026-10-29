function M = chan_model_est(kind, gj, Ts, gin, act, bs, oact, dl)
%CHAN_MODEL_EST Implemented sampled-data channel model of chan_model_inner
%   with a separate measurement delay for each measured signal (used for the
%   estimator-aware covariance). Identical to chan_model_inner except
%   that dl = [d_y d_v d_th d_om] (samples) replaces the common plant delay
%   act(2): the controller receives y(k-d_y), v(k-d_v), theta(k-d_th) and
%   omega(k-d_om) plus noise. dl = act(2)*[1 1 1 1] gives exactly the model of
%   chan_model_inner (checked in cov_estimator_aware_check).
%   act(1) = plant rotor time constant; oact = [tau_m, nd] assumed INSIDE
%   the controller (unchanged, the observer is not retuned).
%   Inputs  w = [w_d, n_y, n_v, n_th, n_om, r, w_u]; outputs y, u, u_t and,
%   for the cascade, the inner virtual control (Ci, Di).
if nargin < 4, gin = []; end
if nargin < 6 || isempty(bs), bs = 1; end
if nargin < 7 || isempty(oact), oact = act; end
if nargin < 8 || isempty(dl), dl = round(act(2))*[1 1 1 1]; end
dl = round(dl(:)');
g = 9.81; casc = ~isempty(gin);
alpha = exp(-Ts/act(1));
alo = exp(-Ts/oact(1)); ndo = round(oact(2));
lk = kind; if kind == 7 || kind == 8 || kind == 9, lk = 4; end
ko = kind; if casc && kind == 7, ko = 4; end
if casc && kind==9,ko=8;end          % outer x-y law
no = max(dl(1:2)); ni = max(dl(3:4));   % history lengths (outer, inner)
nx1 = 4 + 2*no + 16; nx = nx1 + casc*(3 + 2*ni + 16); nw = 7;
A = zeros(nx); B = zeros(nx,nw); Cy = zeros(1,nx); Cu = zeros(1,nx); Ct = zeros(1,nx);
Du = zeros(1,nw); Dt = zeros(1,nw); Ci = zeros(1,nx); Di = zeros(1,nw);
for c = 1:nx, e = zeros(nx,1); e(c) = 1; [A(:,c), Cy(c), Cu(c), Ct(c), Ci(c)] = f1(e, zeros(nw,1)); end
for c = 1:nw, e = zeros(nw,1); e(c) = 1; [B(:,c), ~, Du(c), Dt(c), Di(c)] = f1(zeros(nx,1), e); end
keep = true(1,nx); keep(3) = false;          % d (random walk) is handled separately
for i = 4:nx
    ei = zeros(1,nx); ei(i) = 1;
    if isequal(A(i,:), ei) && isequal(A(:,i)', ei) && ~any(B(i,:)), keep(i) = false; end
end
M.A = A; M.B = B; M.Cy = Cy; M.Cu = Cu; M.Ct = Ct; M.Du = Du; M.Dt = Dt;
M.Ci = Ci; M.Di = Di; M.keep = find(keep); M.id = 3; M.dl = dl;
    function [Xn, yo, uo, to, iuo] = f1(X, w)
        iuo = 0; y = X(1); v = X(2); d = X(3); a = X(4);
        yb = X(5:4+no); vb = X(5+no:4+2*no); cs = X(5+2*no:nx1);
        ym = pick(yb, y, dl(1)) + w(2); vm = pick(vb, v, dl(2)) + w(3);
        u  = chan_output(lk, gj, cs(1:4), ym, vm, w(6), 0, 0, 0);
        ut = u + w(7);
        if ~casc
            csn = chan_update(kind, gj, cs, ym, vm, w(6), u, Ts, inf, alo, ndo);
            yn = y + Ts*v + Ts^2/2*(d + a);  vn = v + Ts*(d + a);
            an = alpha*a + (1-alpha)*bs*ut;  Xi = [];
        else
            o = nx1;
            th = X(o+1); om = X(o+2); ath = X(o+3);
            thb = X(o+4:o+3+ni); omb = X(o+4+ni:o+3+2*ni); ci = X(o+4+2*ni:nx);
            thm = pick(thb, th, dl(3)) + w(4); omm = pick(omb, om, dl(4)) + w(5);
            if ko == 4 || ko == 8, uae = g*thm; else, uae = u; end
            csn = chan_update(ko, gj, cs, ym, vm, w(6), uae, Ts, inf, alo, ndo);
            thd = ut/g;
            ui  = chan_output(lk, gin, ci(1:4), thm, omm, thd, 0, 0, 0);
            iuo = ui; cin = chan_update(kind, gin, ci, thm, omm, thd, ui, Ts, inf, alo, ndo);
            yn = y + Ts*v + Ts^2/2*(d + g*th);  vn = v + Ts*(d + g*th);  an = a;
            thn = th + Ts*om + Ts^2/2*ath;  omn = om + Ts*ath;  athn = alpha*ath + (1-alpha)*bs*ui;
            Xi = [thn; omn; athn; shift(thb, th); shift(omb, om); cin];
        end
        Xn = [yn; vn; d + w(1); an; shift(yb, y); shift(vb, v); csn; Xi];
        yo = y; uo = u; to = ut;
    end
end
function s = pick(b, now, d)
if d >= 1, s = b(d); else, s = now; end
end
function b = shift(b, newest)
if isempty(b), return, end
b = [newest; b(1:end-1)];
end
