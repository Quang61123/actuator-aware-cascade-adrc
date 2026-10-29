function [os, ts] = step_response_metrics(M, Ts, stride, Tend)
%STEP_RESPONSE_METRICS Overshoot [%] and 2% settling time [s] of the unit reference step
%   of a sampled closed-loop channel model M (chan_model or chan_model_est;
%   for the horizontal class the model contains the inner loop, i.e. the
%   cascade). Response requirement: horizon 400 s.
%   stride > 1 evaluates the exact response at every stride-th sample
%   (lifted system, step input held constant; used only to screen during
%   searches; final checks use stride = 1).
if nargin < 3 || isempty(stride), stride = 1; end
if nargin < 4 || isempty(Tend), Tend = 400; end
k = M.keep; A = M.A(k,k); B = M.B(k,6); c = M.Cy(k);
if max(abs(eig(A))) >= 1, os = inf; ts = inf; return, end
if stride > 1
    As = A^stride; Bs = zeros(size(B)); term = B;
    for j = 1:stride, Bs = Bs + term; term = A*term; end
    A = As; B = Bs;
end
h = Ts*stride; N = ceil(Tend/h); t = (0:N-1)'*h;
y = ltitr(A, B, ones(N,1))*c';
os = max(0, 100*(max(y) - 1));
out = find(abs(y - 1) > 0.02, 1, 'last');
if isempty(out), ts = 0; elseif out == N, ts = inf; else, ts = t(out+1); end
end
