function z = cov_estimator_aware(M, NF)
%COV_ESTIMATOR_AWARE Stationary covariances of a channel model driven by
%   coloured (estimator-shaped) measurement noise.
%   M  : model from chan_model_est / chan_model_inner / chan_model
%        (noise inputs are the columns 2..1+p of M.B, p = 4 for a cascade
%        [n_y n_v n_th n_om], p = 2 for a single channel [n_y n_v]).
%   NF : shaping filter of the noise vector n (p outputs), white v ~ N(0,I):
%           z(k+1) = Af z(k) + Bf v(k),   n(k) = Cf z(k) + Df v(k)
%        (Af may be empty: white noise n = Df v with covariance Df*Df').
%   The filter is appended to the channel state:
%     Xa = [X; z],  Fa = [F  Bn*Cf; 0 Af],  Ga = [Bn*Df; Bf],
%     output row ca = [c  dn*Cf], feed-through da = dn*Df,
%   and the white-noise formulas of chan_var / cascade_stats are applied
%   unchanged with unit noise intensity.
%   z = [sigma_e, sigma_u, sigma_du, sigma_du_inner, rho]
%   (tracking, outer virtual control, its increment, inner virtual-control
%   increment of the cascade; rho = spectral radius of the closed loop).
k = M.keep; F = M.A(k,k); p = size(NF.Df,1);
Bn = M.B(k,2:1+p);
nz = size(NF.Af,1); n = numel(k);
Fa = [F, Bn*NF.Cf; zeros(nz,n), NF.Af];
Ga = [Bn*NF.Df; NF.Bf];
rho = max(abs(eig(Fa)));
z = inf(1,5); z(5) = rho;
if rho >= 1, return, end
Pa = dlyap(Fa, Ga*Ga');
[ve] = outvar(M.Cy(k), zeros(1,p));                 % tracking (noise part)
[vu, vdu] = outvar(M.Cu(k), M.Du(2:1+p));
vdui = NaN;
if isfield(M,'Ci') && any(M.Ci(k)), [~, vdui] = outvar(M.Ci(k), M.Di(2:1+p)); end
z = [sqrt(max(0,ve)), sqrt(max(0,vu)), sqrt(max(0,vdu)), sqrt(max(0,vdui)), rho];
    function [v, vd] = outvar(c, dn)
        ca = [c, dn*NF.Cf]; da = dn*NF.Df;
        v = ca*Pa*ca' + da*da';
        % du_k = ca(Fa-I)Xa_{k-1} + (ca*Ga - da) v_{k-1} + da v_k
        h = ca*(Fa - eye(size(Fa,1))); q = ca*Ga - da;
        vd = h*Pa*h' + q*q' + da*da';
    end
end
