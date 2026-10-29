function z = est_stats(co, S, NF)
%EST_STATS Estimator-aware cascade statistics of a design.
%   As cascade_stats, but with the estimator-coloured noise models
%   (pooled over designs) and the signal delays of fit_noise_per_design:
%   z = [sigma_e, sigma_u, sigma_du, sigma_du_inner, rho], where sigma_e
%   adds the random-walk-disturbance part of cascade_stats to the noise
%   part, and sigma_du_inner is the larger of the pitch (x) and roll (y)
%   cascade increments.
ds = NF.dsig; oa = [co.tauhat co.ndhat];
MX = chan_model_est(co.kind, co.G(1:6)', S.Ts, co.G(19:24)', S.act, [], oa, ds([1 4 8 11]));
MY = chan_model_est(co.kind, co.G(1:6)', S.Ts, co.G(19:24)', S.act, [], oa, ds([2 5 7 10]));
zX = cov_estimator_aware(MX, NF.X_pooled); zY = cov_estimator_aware(MY, NF.Y_pooled);
z = inf(1,5); z(5) = max(zX(5), zY(5));
if z(5) >= 1, return, end
k = MX.keep; F = MX.A(k,k); bd = MX.A(k,MX.id); ce = MX.Cy(k); xss = (eye(numel(k)) - F)\bd;
if abs(ce*xss) > 1e-7*max(1, norm(ce)*norm(xss)), return, end
vd = S.qd(1)*S.Ts*(xss'*dlyap(F', ce'*ce)*xss);
z = [sqrt(zX(1)^2 + vd), max(zX(2),zY(2)), max(zX(3),zY(3)), max(zX(4),zY(4)), z(5)];
end
