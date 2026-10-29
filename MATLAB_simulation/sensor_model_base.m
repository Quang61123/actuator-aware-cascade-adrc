function meas = sensor_model_base(x, n, sig)
%SENSOR_MODEL_BASE Noisy measurements. n = 12x1 standard-normal sample,
%   sig = [sp sv sa sg] std of position, velocity, angle, gyro.
%   Euler-angle rates are obtained from the noisy gyro through W^-1(eta_meas).
%#codegen
pos = x(1:3) + sig(1)*n(1:3);
vel = x(4:6) + sig(2)*n(4:6);
eta = x(7:9) + sig(3)*n(7:9);
gyr = x(10:12) + sig(4)*n(10:12);
phi = eta(1); th = eta(2);
Wi = [1, sin(phi)*tan(th), cos(phi)*tan(th); 0, cos(phi), -sin(phi); 0, sin(phi)/cos(th), cos(phi)/cos(th)];
meas = [pos; vel; eta; Wi*gyr];
end


