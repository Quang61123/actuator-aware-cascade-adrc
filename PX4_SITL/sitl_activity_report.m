function S = sitl_activity_report()
%SITL_ACTIVITY_REPORT Comparison of PX4 SITL logs with the MATLAB results.
%   Reads ../results/sitl_activity/runs/*.mat (run_sitl_activity) and writes to
%   ../results/sitl_activity/. Two controller input modes:
%     EKF  (tags act_*)    : controller fed by EKF2 and the filtered gyro (nominal setup)
%     raw  (tags act_raw_*) : estimator bypass, controller fed by the
%                           measurement model (nominal white noise, 5-ms delay)
%   Files (prefix act_ for EKF, act_raw_ for raw):
%     *_s6_per_run.csv, *_hover_per_run.csv  per-run values
%     *_s6_summary.csv    S6, seeds 91001-91010: RMSE and command increment
%                         J on the unquantized command (J_unq,
%                         PX4 float output) and on the quantized command
%                         received by the plant (J_q), next to the MATLAB
%                         E1/E0 increment of the same seeds (trials.csv)
%     *_hover_summary.csv noise-only hover, seeds 96001-96010, window 20-60 s:
%                         sigma_Delta u (roll, pitch; pooled as in hover_per_seed.csv),
%                         next to the MATLAB measurement (hover_per_seed.csv)
%                         and the linear prediction
%     act_estimator.csv    EKF2 estimate error (RMS, SD, delay) per design/case
%     act_psd.png          mean PSD of the x-position and roll estimate errors (EKF hover)
%     act_ordering.txt     design ordering of each activity metric (both modes)
%     act_compare.csv      MATLAB vs SITL-EKF vs SITL-raw, one row per design
%     act_raw_measurement.csv  noise SD, increment SD and lag of the measurement
%                         seen by the controller in raw mode (check of the bypass)
%   MATLAB values are only read; nothing in results is recomputed.
here = fileparts(mfilename('fullpath'));
RP = fullfile(here,'..','results'); RACT = fullfile(RP,'sitl_activity'); RR = fullfile(RACT,'runs');
names = {1,'PID-Budget'; 2,'AA-Budget'; 4,'AA-HalfGain'; 8,'PID-Budget-Strict'; 9,'AA-Budget-Strict'; 10,'AA-HalfGain-Strict'};
short = containers.Map({'PID-Budget','AA-Budget','AA-HalfGain','PID-Budget-Strict','AA-Budget-Strict','AA-HalfGain-Strict'}, ...
                       {'PID-Budget','AA-Budget','AA-HalfGain','PID-Budget-S','AA-Budget-S','AA-HalfGain-S'});
pred = containers.Map({1,2,4,8,9,10}, {1.380,0.937,1.023,0.592,0.496,0.521});   % predicted cascade sigma_Delta u (cascade_stats)
TR = readtable(fullfile(RP,'trials.csv'));
HV = readtable(fullfile(RP,'hover_per_seed.csv'));
f = dir(fullfile(RR,'act*.mat'));
A = table(); PSDp = {}; PSDa = {}; MZ = table();
fq = 0.25:0.25:45;                                         % common grid below the 50-Hz Nyquist of the 100-Hz position log
for i = 1:numel(f)
    D = load(fullfile(RR,f(i).name)); r = D.run; X = D.X; M = D.M;
    raw = isfield(r,'raw') && r.raw;
    nm = names{[names{:,1}] == r.ctrl, 2};
    if r.hover, cs = "hover"; else, cs = "S6"; end
    if raw, md = "raw"; else, md = "EKF"; end
    row = table(md, cs, string(nm), r.seed, M.rmse3*1e3, X.J_q, X.J_unq, X.frac_q_nonzero, ...
        X.sd_du(1), X.sd_du(2), X.sd_du(3), X.sd_du(4), X.sd_u(1), X.sd_u(2), ...
        X.est_pos_rms(1), X.est_pos_rms(2), X.est_pos_rms(3), X.est_vel_rms(1), X.est_vel_rms(2), X.est_vel_rms(3), ...
        X.est_att_rms(1), X.est_att_rms(2), X.est_att_rms(3), X.est_att_sd(1), X.est_att_sd(2), X.est_pos_sd(1), X.est_rate_hz, X.lag_pos_ms, X.lag_att_ms, ...
        X.n_motor_samples, X.n_expected, X.frac_consecutive, X.max_gap_ms, r.miss, double(r.diverged), ...
        'VariableNames', {'Mode','Case','Ctrl','Seed','RMSE_mm','J_q','J_unq','FracQuantNonzero', ...
        'SdDu_roll','SdDu_pitch','SdDu_yaw','SdDu_z','SdU_roll','SdU_pitch', ...
        'EstPosRMS_x','EstPosRMS_y','EstPosRMS_z','EstVelRMS_x','EstVelRMS_y','EstVelRMS_z', ...
        'EstAttRMS_roll','EstAttRMS_pitch','EstAttRMS_yaw','EstAttSD_roll','EstAttSD_pitch','EstPosSD_x','EstRate_Hz','EstLagPos_ms','EstLagAtt_ms', ...
        'MotorSamples','MotorExpected','FracConsecutive','MaxGap_ms','Miss','Diverged'});
    A = [A; row]; %#ok<AGROW>
    PSDp{end+1,1} = interp1(X.psd_f(:), X.psd_pos(:), fq); PSDa{end+1,1} = interp1(X.psd_f(:), X.psd_att(:), fq); %#ok<AGROW>
    if raw && isfield(X,'meas_err_sd_x')
        MZ = [MZ; table(cs, string(nm), r.seed, X.meas_err_sd_x, X.meas_err_sd_phi, X.meas_inc_sd(1), X.meas_inc_sd(2), ...
            X.meas_lag_ms(1), X.meas_lag_ms(2), 'VariableNames', {'Case','Ctrl','Seed','MeasErrSD_x','MeasErrSD_roll', ...
            'MeasIncSD_x','MeasIncSD_roll','MeasLag_x_ms','MeasLag_roll_ms'})]; %#ok<AGROW>
    end
end
assert(~isempty(A), 'sitl_activity_report: no runs found');
if ~isempty(MZ), writetable(MZ, fullfile(RACT,'act_raw_measurement.csv')); end
modes = ["EKF","raw"]; pre = ["act","act_raw"]; S6m = struct(); HSm = struct();
fid = fopen(fullfile(RACT,'act_ordering.txt'),'w');
ord = @(T,v) strjoin(T.Design(sortrows(table((1:height(T))', T.(v)),2).Var1), ' < ');
primSet = ["PID-Budget","AA-Budget","AA-HalfGain"];
for m = 1:2
    B = A(A.Mode == modes(m),:);
    if isempty(B), continue; end
    writetable(B(B.Case=="S6",:), fullfile(RACT,pre(m) + "_s6_per_run.csv"));
    writetable(B(B.Case=="hover",:), fullfile(RACT,pre(m) + "_hover_per_run.csv"));
    % ---------------- S6 summary ----------------
    S6 = table();
    for j = 1:size(names,1)
        c = names{j,1}; nm = names{j,2}; a = B(B.Case=="S6" & B.Ctrl==nm,:); a = sortrows(a,'Seed');
        if isempty(a), continue; end
        ex = 'E1'; if c >= 8, ex = 'E0'; end
        t = TR(strcmp(TR.Experiment,ex) & strcmp(TR.Ctrl,nm) & ismember(TR.Seed,a.Seed),:); t = sortrows(t,'Seed');
        mI = NaN; mR = NaN; if height(t) == height(a), mI = mean(t.Increment); mR = mean(t.RMSE_mm); end
        S6 = [S6; table(string(short(nm)), height(a), mR, mean(a.RMSE_mm), std(a.RMSE_mm), mI, mean(a.J_unq), std(a.J_unq), ...
            mean(a.J_q), std(a.J_q), mean(a.FracQuantNonzero), mean(a.SdDu_roll), mean(a.SdDu_pitch), sum(a.Miss > 0), sum(a.Diverged), ...
            'VariableNames', {'Design','N','MATLAB_RMSE_mm','SITL_RMSE_mm','SITL_RMSE_sd','MATLAB_J','SITL_J_unq','SITL_J_unq_sd', ...
            'SITL_J_q','SITL_J_q_sd','FracQuantNonzero','SITL_SdDu_roll','SITL_SdDu_pitch','RunsWithMiss','Diverged'})]; %#ok<AGROW>
    end
    writetable(S6, fullfile(RACT,pre(m) + "_s6_summary.csv"));
    % ---------------- hover summary ----------------
    HS = table();
    for j = 1:size(names,1)
        c = names{j,1}; nm = names{j,2}; a = B(B.Case=="hover" & B.Ctrl==nm,:);
        if isempty(a), continue; end
        h = HV(strcmp(HV.Ctrl,nm),:);
        mr = sqrt(mean(h.MeasuredSigmaDU.^2)); mp = sqrt(mean(h.PitchMeasuredSigmaDU.^2));
        sr = sqrt(mean(a.SdDu_roll.^2)); sp = sqrt(mean(a.SdDu_pitch.^2));
        HS = [HS; table(string(short(nm)), height(a), pred(c), mr, mp, max(mr,mp), sr, sp, max(sr,sp), max(sr,sp)/pred(c), ...
            mean(a.J_unq), mean(a.J_q), mean(a.FracQuantNonzero), mean(a.RMSE_mm), ...
            'VariableNames', {'Design','N','Pred_SigmaDU','MATLAB_roll','MATLAB_pitch','MATLAB_max','SITL_roll','SITL_pitch','SITL_max', ...
            'SITL_over_Pred','SITL_J_unq','SITL_J_q','FracQuantNonzero','SITL_RMSE_mm'})]; %#ok<AGROW>
    end
    writetable(HS, fullfile(RACT,pre(m) + "_hover_summary.csv"));
    S6m.(modes(m)) = S6; HSm.(modes(m)) = HS;
    % ---------------- ordering ----------------
    fprintf(fid,'===== SITL input mode: %s =====\n', modes(m));
    for g = {ismember(S6.Design,primSet), ~ismember(S6.Design,primSet)}
        T = S6(g{1},:); if isempty(T), continue; end
        fprintf(fid,'S6 MATLAB J      : %s\n', ord(T,'MATLAB_J'));
        fprintf(fid,'S6 SITL J_unq    : %s\n', ord(T,'SITL_J_unq'));
        fprintf(fid,'S6 SITL J_q      : %s\n', ord(T,'SITL_J_q'));
        fprintf(fid,'S6 SITL sd_du rl : %s\n\n', ord(T,'SITL_SdDu_roll'));
    end
    for g = {ismember(HS.Design,primSet), ~ismember(HS.Design,primSet)}
        T = HS(g{1},:); if isempty(T), continue; end
        fprintf(fid,'hover pred sigma_du : %s\n', ord(T,'Pred_SigmaDU'));
        fprintf(fid,'hover MATLAB max    : %s\n', ord(T,'MATLAB_max'));
        fprintf(fid,'hover SITL max      : %s\n', ord(T,'SITL_max'));
        fprintf(fid,'hover SITL J_unq    : %s\n\n', ord(T,'SITL_J_unq'));
    end
end
fclose(fid);
% ---------------- three-way comparison ----------------
if isfield(S6m,'EKF') && isfield(S6m,'raw') && isfield(HSm,'EKF') && isfield(HSm,'raw')
    a = S6m.EKF; b = S6m.raw; h = HSm.EKF; g = HSm.raw;
    [~,ib] = ismember(a.Design, b.Design); [~,ih] = ismember(a.Design, h.Design); [~,ig] = ismember(a.Design, g.Design);
    C = table(a.Design, a.MATLAB_RMSE_mm, a.SITL_RMSE_mm, b.SITL_RMSE_mm(ib), a.MATLAB_J, a.SITL_J_unq, b.SITL_J_unq(ib), b.SITL_J_q(ib), ...
        h.Pred_SigmaDU(ih), h.MATLAB_max(ih), h.SITL_max(ih), g.SITL_max(ig), g.SITL_over_Pred(ig), ...
        'VariableNames', {'Design','MATLAB_RMSE_mm','EKF_RMSE_mm','Raw_RMSE_mm','MATLAB_J','EKF_J_unq','Raw_J_unq','Raw_J_q', ...
        'Pred_SigmaDU','MATLAB_hover_SigmaDU','EKF_hover_SigmaDU','Raw_hover_SigmaDU','Raw_over_Pred'});
    writetable(C, fullfile(RACT,'act_compare.csv'));
else
    C = table();
end
% ---------------- estimator ----------------
E = groupsummary(A, {'Mode','Case','Ctrl'}, 'mean', {'EstPosRMS_x','EstPosRMS_y','EstPosRMS_z','EstVelRMS_x','EstVelRMS_y', ...
    'EstAttRMS_roll','EstAttRMS_pitch','EstAttRMS_yaw','EstAttSD_roll','EstAttSD_pitch','EstPosSD_x','EstRate_Hz','EstLagPos_ms','EstLagAtt_ms','MaxGap_ms'});
writetable(E, fullfile(RACT,'act_estimator.csv'));
% ---------------- PSD figure (EKF hover runs) ----------------
ih = find(A.Case=="hover" & A.Mode=="EKF");
if ~isempty(ih)
    fig = figure('Visible','off','Position',[100 100 800 320]);
    subplot(1,2,1); loglog(fq, mean(cell2mat(PSDp(ih)),1)); grid on
    xlabel('f [Hz]'); ylabel('PSD x-est. error [m^2/Hz]'); title('EKF2 position error (hover)');
    subplot(1,2,2); loglog(fq, mean(cell2mat(PSDa(ih)),1)); grid on
    xlabel('f [Hz]'); ylabel('PSD roll-est. error [rad^2/Hz]'); title('EKF2 roll error (hover)');
    exportgraphics(fig, fullfile(RACT,'act_psd.png'), 'Resolution', 150); close(fig);
end
S = struct('S6', S6m, 'hover', HSm, 'compare', C, 'est', E, 'meas', MZ, 'all', A);
end
