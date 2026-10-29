function [os, ts, itae] = chan_step(kind, gj, Ts, gin, act, Tend, oact)
%CHAN_STEP Reference-step overshoot [%], 2% settling time [s] and an
%   ITAE cost (reference step + unit disturbance step) of the channel model.
if nargin < 7, oact = []; end
M = chan_model(kind, gj, Ts, gin, act, [], oact);
k = M.keep; A = M.A(k,k); c = M.Cy(k);
if max(abs(eig(A))) >= 1, os = inf; ts = inf; itae = inf; return, end
N = round(Tend/Ts); t = (0:N-1)'*Ts;
y = ltitr(A, M.B(k,6), ones(N,1))*c';
os = max(0, 100*(max(y) - 1));
out = find(abs(y - 1) > 0.02, 1, 'last');
if isempty(out), ts = 0; elseif out == N, ts = inf; else, ts = t(out+1); end
if nargout > 2
    yd = ltitr(A, M.A(k,M.id), ones(N,1))*c';       % response to a unit d step
    itae = Ts*sum(t.*abs(1 - y)) + Ts*sum(t.*abs(yd));
end
end


