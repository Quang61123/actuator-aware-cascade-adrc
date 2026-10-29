function T = sitl_activity_sensors()
%SITL_ACTIVITY_SENSORS Noise seen by the controller in PX4 SITL.
%   For each noise-only hover ulog (ulog_act/act_hover_*.ulg) the 1-ms sample
%   increments of the signals the controller reads are compared with the raw
%   sensor signal. During hover the true motion changes little within 1 ms,
%   so the increment standard deviation is dominated by noise; for white
%   noise of standard deviation sigma it equals sqrt(2)*sigma, i.e. the value
%   seen by the MATLAB controller.
%     gyro : vehicle_angular_velocity (PX4 filtered, IMU_GYRO_CUTOFF) vs
%            sensor_combined.gyro_rad (raw simulated gyro), body x/y axes
%     att  : vehicle_attitude roll/pitch (EKF2) increments
%     pos  : vehicle_local_position x/y increments (logged at 100 Hz, so the
%            increments are over 10 ms; reported for completeness)
%   Window: 20-60 s after the first armed actuator_motors sample.
%   Output: ../results/sitl_activity/act_controller_input_noise.csv
here = fileparts(mfilename('fullpath'));
RACT = fullfile(here,'..','results','sitl_activity'); UW = fullfile(px4sitl_work_dir(),'ulog_act');
f = dir(fullfile(UW,'act_hover_*.ulg'));
T = table();
for i = 1:numel(f)
    u = ulogreader(fullfile(UW,f(i).name));
    d = readTopicMsgs(u,'TopicNames',{'actuator_motors','vehicle_angular_velocity','sensor_combined','vehicle_attitude','vehicle_local_position'});
    nm = string(d.TopicNames); g = @(s) d.TopicMessages{nm == s};
    A = g("actuator_motors"); c = double(A.control); ta = seconds(A.timestamp);
    t0 = ta(find(c(:,10) > 0, 1));                       % first armed controller sample
    win = t0 + [20 60];
    G = g("vehicle_angular_velocity"); [sdF, rF] = incsd(seconds(G.timestamp), double(G.xyz(:,1:2)), win, 1e-3);
    S = g("sensor_combined");          [sdR, rR] = incsd(seconds(S.timestamp), double(S.gyro_rad(:,1:2)), win, 1e-3);
    Q = g("vehicle_attitude"); q = double(Q.q);
    eul = [atan2(2*(q(:,1).*q(:,2)+q(:,3).*q(:,4)), 1-2*(q(:,2).^2+q(:,3).^2)), asin(max(-1,min(1,2*(q(:,1).*q(:,3)-q(:,4).*q(:,2)))))];
    [sdA, rA] = incsd(seconds(Q.timestamp), eul, win, 1e-3);
    P = g("vehicle_local_position"); [sdP, rP] = incsd(seconds(P.timestamp), [double(P.x) double(P.y)], win, 1e-2);
    tok = regexp(f(i).name,'act_hover_c(\d+)_seed(\d+)','tokens','once');
    T = [T; table(str2double(tok{1}), str2double(tok{2}), mean(sdR), mean(sdF), mean(sdF)/mean(sdR), mean(sdA), mean(sdP), ...
        rR, rF, rA, rP, 'VariableNames', {'Ctrl','Seed','GyroRawIncSD','GyroFiltIncSD','GyroFiltOverRaw', ...
        'AttIncSD_1ms','PosIncSD_10ms','FracRaw','FracFilt','FracAtt','FracPos'})]; %#ok<AGROW>
end
T = sortrows(T, {'Ctrl','Seed'});
writetable(T, fullfile(RACT,'act_controller_input_noise.csv'));
end

function [sd, frac] = incsd(t, y, win, dt)
% standard deviation of increments over consecutive samples spaced dt
k = t >= win(1) & t < win(2); t = t(k); y = y(k,:);
[t, iu] = unique(t); y = y(iu,:);
ok = abs(diff(t) - dt) < 0.25*dt; frac = mean(ok);
dy = diff(y); sd = std(dy(ok,:), 1, 1);
end
