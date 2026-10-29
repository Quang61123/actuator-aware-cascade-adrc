function m = run_metrics(L, win, pv)
%RUN_METRICS Window metrics of one simulation run (Inf if the run diverged).
k = L.t >= win(1) & L.t < win(2);
bad = any(~isfinite(L.x(:))) || max(abs(L.x(:,1:3)),[],'all') > 50;
e = L.r(k,1:3) - L.x(k,1:3); en = sqrt(sum(e.^2,2));
m.rmse3 = sqrt(mean(en.^2)); m.max3 = max(en);
m.rmsez = sqrt(mean(e(:,3).^2)); m.rmsepsi = sqrt(mean((L.r(k,4)-L.x(k,9)).^2));
dwc = diff(L.wcmd(k,:)); m.chat_cmd = sqrt(mean(dwc(:).^2));      % command increments [rad/s]
dw  = diff(L.x(k,13:16)); m.chat_act = sqrt(mean(dw(:).^2));      % actual rotor-speed increments
wh = sqrt(pv(1)*pv(2)/(4*pv(4)));
m.effort = sqrt(mean((L.x(k,13:16) - wh).^2, 'all'));
% rotor-speed fluctuation about its window mean (oscillation indicator, threshold 20 rad/s)
m.rotor_std = sqrt(mean(var(L.x(k,13:16), 1, 1)));
% saturation duty cycles (fraction of samples in the window)
wmax = pv(11); tol = 1e-6;
m.sat_rotor = mean(any(L.wcmd(k,:) >= wmax - tol | L.wcmd(k,:) <= tol, 2));
if isfield(L,'Tc')
    angc = any(abs(L.aref(k,1:2)) >= 0.35 - 1e-9, 2);
    Tcl  = L.Tc(k) <= tol | L.Tc(k) >= pv(12) - tol;
    tcl  = any(abs(L.tauc(k,:)) >= pv(13:15)' - tol, 2);
    m.sat_virt = mean(angc | Tcl | tcl);
else, m.sat_virt = NaN; end
if bad, f = fieldnames(m); for q = 1:numel(f), m.(f{q}) = Inf; end, end
m.diverged = bad;
end


