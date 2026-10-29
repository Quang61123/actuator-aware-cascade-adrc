function pv = quad_params(scale, tau_m, nd)
%QUAD_PARAMS Base parameter vector + rotor time constant and loop delay.
%   pv(1:15) as quad_params_base; pv(16) = tau_m [s]; pv(17) = nd [samples].
%   scale = [s_m s_I s_A].  Defaults: tau_m = 0.03 s, nd = 5.
%   tau_m and nd are modelling ASSUMPTIONS (not given in the
%   baseline references); they are swept in sweep and robust-grid experiments.
if nargin < 1 || isempty(scale), scale = [1 1 1]; end
if nargin < 2 || isempty(tau_m), tau_m = 0.03; end
if nargin < 3 || isempty(nd), nd = 5; end
if ~(isscalar(tau_m) && isfinite(tau_m) && tau_m > 0), error('quad_params:tau','tau_m must be a positive scalar'); end
if ~(isscalar(nd) && nd == round(nd) && nd >= 0 && nd <= 10)
    error('quad_params:nd','nd must be an integer in 0..10 (sensor buffer size)');
end
pv = [quad_params_base(scale); tau_m; nd];
end



