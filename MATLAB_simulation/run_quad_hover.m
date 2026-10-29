function L = run_quad_hover(kind, G, scen, pv_true, pv_nom, Tend, Ts, E, sig, ndSeries)
%RUN_QUAD_HOVER Additional stress protocol; shared RK4 plant/controller: RK4 plant (16 states incl. rotor
%   dynamics), controller at Ts, measurement delay of pv_true(17) samples,
%   exogenous noise / gyro bias / wind from make_exo (struct E).
nsub=1;
assert(Ts>0 && abs(Tend/Ts-round(Tend/Ts))<1e-8,'Invalid simulation horizon');
N = round(Tend/Ts);ndSeries=ndSeries(:);assert(numel(ndSeries)==N&&all(ndSeries==round(ndSeries))&&all(ndSeries>=0&ndSeries<=30),'Delay must be an integer sequence in 0..30');
assert(isequal(size(E.noise),[12 N]) && isequal(size(E.bias),[3 N]) && isequal(size(E.wind),[3 N]),'Exogenous signal size mismatch');
x = zeros(16,1); x(13:16) = sqrt(pv_true(1)*pv_true(2)/(4*pv_true(4)));  % hover rotor speed
cs = zeros(96,1); buf = zeros(12, 31); init = 0;   % measurement history (col 1 = newest)
L.t = (0:N-1)'*Ts; L.x = zeros(N,16); L.r = zeros(N,10); L.aref = zeros(N,3);
L.T = zeros(N,1); L.tau = zeros(N,3); L.wcmd = zeros(N,4);
L.u=zeros(N,6); L.Tc = zeros(N,1); L.tauc = zeros(N,3);      % controller outputs before allocation 
for i = 1:N
    t = (i-1)*Ts;
    r=zeros(10,1);
    mnow=sensor_model(x,E.noise(:,i),sig,E.bias(:,i));
    if init==0,buf=repmat(mnow,1,31);end
    buf=[mnow,buf(:,1:30)];init=1;meas=buf(:,ndSeries(i)+1);
    [T, tau, aref, cs, uv] = ctrl_law_logged(kind, G, cs, meas, r, pv_nom, Ts);
    [wc, Ta, taua] = quad_mixer(T, tau, pv_true);
    Fd=zeros(3,1);Td=zeros(3,1);
    L.u(i,:)=uv';L.x(i,:) = x'; L.r(i,:) = r'; L.aref(i,:) = aref'; L.T(i) = Ta; L.tau(i,:) = taua'; L.wcmd(i,:) = wc'; L.Tc(i) = T; L.tauc(i,:) = tau';
    for sub=1:nsub
    h=Ts/nsub;
    k1 = quad_dynamics(x,         wc, Fd, Td, pv_true);
    k2 = quad_dynamics(x+h/2*k1, wc, Fd, Td, pv_true);
    k3 = quad_dynamics(x+h/2*k2, wc, Fd, Td, pv_true);
    k4 = quad_dynamics(x+h*k3,   wc, Fd, Td, pv_true);
    x = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    if any(~isfinite(x))||max(abs(x(1:3)))>50
      L.x(i+1:end,:)=NaN;L.wcmd(i+1:end,:)=NaN;break
    end
end
end



