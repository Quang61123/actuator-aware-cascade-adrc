function R = sitl_log_extract(ulogFile)
%SITL_LOG_EXTRACT Controller-rate signals of one SITL run.
%   Reads debug_array from the PX4 ulog (written by the code-generated
%   controller at every 1-ms step, see sitl_ctrl_step). Each message holds
%   the block of the current and of the previous step, so a message dropped
%   by the logger is recovered from the next one. Output:
%     R.tc   controller time k*Ts [s] (k = 0 at the first armed step)
%     R.meas 12 measured signals used by ctrl_law (model frame:
%            x y z, vx vy vz, phi theta psi, Euler-angle rates)
%     R.ua   [u_phi u_theta u_psi u_z] (limited virtual controls)
%     R.wn   unquantized w_cmd/w_max (4 rotors), R.T thrust command
%     R.src  input source (0 EKF2, 1 raw measurement)
%   Only armed samples are kept, duplicates removed and sorted by tc.
%   R.ok marks samples whose predecessor is the previous controller step.
u = ulogreader(ulogFile);
d = readTopicMsgs(u,'TopicNames',{'debug_array'});
D = double(d.TopicMessages{1}.data);
D = D(D(:,1) > 0,:);                             % id > 0: armed samples
B = [D(:,3:24); D(:,25:46)];                     % current and previous blocks
B = B(B(:,1) > 0 | any(B(:,2:13) ~= 0, 2), :);   % drop the empty block of the first step
k = round(B(:,1)*1000);
[k, iu] = unique(k); B = B(iu,:);
R.tc = k/1000; R.meas = single(B(:,2:13)); R.ua = single(B(:,14:17)); R.wn = single(B(:,18:21));
R.T = single(B(:,22)); R.src = median(D(:,2)); R.id = median(D(:,1));
R.ok = [false; diff(k) == 1];
R.n = numel(k); R.frac_ok = mean(R.ok);
end
