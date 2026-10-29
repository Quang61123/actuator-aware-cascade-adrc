function meas = sensor_model(x, n, sig, bias)
%SENSOR_MODEL Sensor model (white noise) + gyro bias (12x1 measurement).
%#codegen
x12 = x(1:12); x12(10:12) = x12(10:12) + bias(:);
meas = sensor_model_base(x12, n, sig);
end


