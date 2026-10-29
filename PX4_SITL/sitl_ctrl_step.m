function [u, dbg] = sitl_ctrl_step(pos, vel, q, gyro, armed, ctrlId, scen, GT, KT, TH, opos, ovel, oq, orate, src)
%SITL_CTRL_STEP One 1-ms step of the cascade controller inside PX4.
%   Runs in the px4_simulink_app module (code generated from the Simulink
%   model built by create_sitl_controller.m). Inputs are the PX4 estimates:
%     pos, vel : vehicle_local_position x,y,z / vx,vy,vz   (NED, m, m/s)
%     q        : vehicle_attitude.q  (FRD body -> NED, [w x y z])
%     gyro     : vehicle_angular_velocity.xyz (FRD body, rad/s)
%     armed    : actuator_armed.armed
%     ctrlId   : controller index (1 PID-Budget, 2 AA-Budget, 3 RO-HalfGain,
%                4 AA-HalfGain, 5-7 A1-A3, 8-10 strict designs),
%                PX4 parameter MC_ROLL_P
%     scen     : scenario 1..6 (reference generator), PX4 parameter MC_PITCH_P
%     GT,KT,TH : 37xNc gain table, 1xNc kinds, 2xNc [tauhat; ndhat]
%   Output u (1x12 single) = actuator_motors.control, u_i = w_cmd,i / w_max
%   (i = 1..4, float, i.e. before the 1/1000 output quantization of the
%   PX4 PWM/simulation output). Logging only (no change of gains or
%   control logic): channels 5..12 carry diagnostics that are logged with
%   actuator_motors in the ulog and are not mapped to any output:
%   u(5:7) = tau./I (limited inner virtual controls u_phi, u_theta, u_psi),
%   u(8) = T*cos(phi)cos(theta)/m - g (limited altitude virtual input),
%   u(9) = controller time k*Ts [s], u(10) = controller id,
%   u(11) = measured x and u(12) = measured roll angle as used by the
%   controller (model frame).
%   Logging only: dbg (1x58 single) = debug_array.data, published
%   at every controller step: dbg(1) = id, dbg(2) = src, dbg(3:24) = block
%   of the current step, dbg(25:46) = block of the previous step (so that a
%   message dropped by the logger is recovered from the next one). Block =
%   [k*Ts, meas(1:12), u(5:8), u(1:4), T]: meas is the measurement vector
%   used by ctrl_law (model frame: x y z, vx vy vz, phi theta psi,
%   Euler-angle rates), u(5:8) the limited virtual controls, u(1:4) the
%   unquantized w_cmd/w_max and T the thrust command.
%   scen = 0: hover at the origin (zero reference), noise-only tests.
%   Diagnostic estimator bypass ( src = PX4 parameter MC_YAWRATE_FF
%   > 0.5, unused by the disabled PX4 controllers): the controller reads
%   vehicle_mocap_odometry instead of the EKF2 and gyro-filter outputs; the
%   simulator then sends every 1 ms the output of sensor_model
%   (nominal white noise on position, velocity, Euler angles and gyro
%   rate, gyro bias, 5-ms delay):
%     opos, ovel : position / velocity (NED), oq : attitude quaternion,
%     orate      : body rate (FRD).
%   The model frame is x north, y west, z up (FLU body); PX4 uses NED/FRD.
%   The control law, mixer and reference are the unchanged functions
%   ctrl_law, quad_mixer and reference_gen (MATLAB_simulation folder).
%#codegen
persistent cs k prv
if isempty(cs), cs = zeros(96,1); k = 0; prv = zeros(1,22,'single'); end
Ts = 1e-3;
u = zeros(1,12,'single'); u(:) = single(NaN);   % NaN = motors stopped (disarmed)
dbg = zeros(1,58,'single');
id = round(double(ctrlId)); if id < 1 || id > size(GT,2), id = 2; end
sc = round(double(scen));   if sc < 0 || sc > 6, sc = 6; end
G = GT(:,id); kind = KT(id);
pv = [quad_params_base([1 1 1]); TH(1,id); TH(2,id)];
if ~armed
    cs(:) = 0; k = 0; prv(:) = 0;
    return
end
pp = double(pos(:)); vv = double(vel(:)); qv = double(q(:)); gv = double(gyro(:));
if double(src(1)) > 0.5   % estimator bypass: raw measurement of the sensor model (sensor_model)
    pp = double(opos(:)); vv = double(ovel(:)); qv = double(oq(:)); gv = double(orate(:));
end
% ---- PX4 (NED/FRD) -> model frame (NWU/FLU) ----
qw = qv(1); qx = qv(2); qy = qv(3); qz = qv(4);
phi_n = atan2(2*(qw*qx + qy*qz), 1 - 2*(qx^2 + qy^2));
th_n  = asin(max(-1, min(1, 2*(qw*qy - qz*qx))));
psi_n = atan2(2*(qw*qz + qx*qy), 1 - 2*(qy^2 + qz^2));
phi = phi_n; th = -th_n; psi = -psi_n;
p = gv(1); qq = -gv(2); r = -gv(3);
sf = sin(phi); cf = cos(phi); tt = tan(th); ct = cos(th);
eta_dot = [1, sf*tt, cf*tt; 0, cf, -sf; 0, sf/ct, cf/ct]*[p; qq; r];
meas = [pp(1); -pp(2); -pp(3); ...
        vv(1); -vv(2); -vv(3); ...
        phi; th; psi; eta_dot];
t = k*Ts;
if sc == 0, rr = zeros(10,1); else, rr = reference_gen(t, sc); end
[T, tau, ~, cs] = ctrl_law(kind, G, cs, meas, rr, pv, Ts);
[w, Ta, ~] = quad_mixer(T, tau, pv);
u(:) = single(0); u(1:4) = single(w(:)'/pv(11));
I3 = pv(7:9); cc = cos(meas(7))*cos(meas(8));
u(5:7) = single(tau(:)'./I3(:)'); u(8) = single(T*cc/pv(1) - pv(2));
u(9) = single(t); u(10) = single(id); u(11) = single(meas(1)); u(12) = single(meas(7));
blk = [single(t), single(meas(:)'), u(5:8), u(1:4), single(T)];
dbg(1) = single(id); dbg(2) = single(src(1)); dbg(3:24) = blk; dbg(25:46) = prv; prv = blk;
k = k + 1;
end
