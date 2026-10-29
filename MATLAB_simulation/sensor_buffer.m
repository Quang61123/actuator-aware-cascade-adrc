function [meas, bufn, initn] = sensor_buffer(x, n, sig, bias, buf, init, nd)
%SENSOR_BUFFER Noisy measurement with a delay of nd samples (nd <= 10).
%   buf (12x11) holds the newest measurement in column 1. At start-up the
%   buffer is filled with the first measurement.
%#codegen
mnow = sensor_model(x, n, sig, bias);
if init == 0, buf = repmat(mnow, 1, 11); end
bufn = [mnow, buf(:, 1:10)];
meas = bufn(:, round(nd) + 1);
initn = 1;
end


