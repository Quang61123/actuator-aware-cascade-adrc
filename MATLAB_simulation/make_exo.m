function E = make_exo(scen, seed, Tend, Ts, S)
%MAKE_EXO Exogenous signals shared by the MATLAB and Simulink runs:
%   E.noise (12xN) standard-normal sensor noise (scaled in sensor_model)
%   E.bias  (3xN)  gyro-bias random walk [rad/s]
%   E.wind  (3xN)  wind force A*v_wind [N] (mean wind + turbulence)
%   Turbulence: per-axis first-order Gauss-Markov (simplified turbulence
%   approximation, time constant L/V), sigma and mean wind in S.
%   The vertical turbulence component is scaled by 0.5 (S.turb_vert_scale),
%   i.e. sigma_w = 0.5 m/s horizontally and 0.25 m/s vertically.
%   The drag coefficient A is taken from quad_params_base (0.25 kg/s).
N = round(Tend/Ts); rng(seed,'twister');   % generator fixed: batch workers default to threefry
E.noise = zeros(12,N); E.bias = zeros(3,N); E.wind = zeros(3,N);
noisy = scen >= 2;
if noisy
    E.noise = randn(12,N);
    b = zeros(3,1); sb = S.gyro_bias_rw*sqrt(Ts);
    for k = 1:N, E.bias(:,k) = b; b = b + sb*randn(3,1); end
end
if scen == 5 || scen == 6
    a = exp(-Ts/S.turb_tc); s = S.turb_sigma; v = zeros(3,1);
    vw = zeros(3,N);
    for k = 1:N
        vw(:,k) = v; v = a*v + s*sqrt(1-a^2)*randn(3,1);
    end
    t = (0:N-1)*Ts; on = double(t >= S.wind_on);
    vmean = S.wind_mean(:)*on;
    vs = 0.5; if isfield(S,'turb_vert_scale'), vs = S.turb_vert_scale; end
    vtot = vmean + vw.*[1;1;vs].*on;            % weaker vertical turbulence
    pvb = quad_params_base(); Adrag = pvb(10);  % nominal linear drag A [kg/s]
    E.wind = Adrag*vtot;                        % wind force entering -A(v - v_w)
end
end
