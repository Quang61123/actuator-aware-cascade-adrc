function mdl = create_sitl_controller()
%CREATE_SITL_CONTROLLER Build the Simulink flight-controller model that is
%   code-generated into PX4 (UAV Toolbox Support Package for PX4, board
%   'PX4 Host Target'). All controllers are compiled into
%   one binary (ids 1-10: 7 nominal/ablation + 3 strict designs; further ids below); the active controller and scenario are selected at run time
%   by the PX4 parameters MC_ROLL_P (controller id) and MC_PITCH_P (scenario),
%   which are unused because the PX4 multicopter controllers are disabled.
%   The parameter MC_YAWRATE_FF (also unused) selects the input
%   source, 0 = EKF2 estimates and filtered gyro (nominal setup), 1 = raw
%   measurement from vehicle_mocap_odometry (estimator-bypass diagnostic).
here = fileparts(mfilename('fullpath'));
simRoot = fullfile(here,'..','MATLAB_simulation');
addpath(here, simRoot);
D = load(fullfile(simRoot,'..','results','design_gains.mat'));   % gains (unchanged)
src = [num2cell(D.C(1:7)), num2cell(D.Strict(1:3))];
% ids 11-12: initial equal-activity designs for the estimator-fed loop, ids 11-12
fE = fullfile(simRoot,'..','results','estimator_noise','est_designs_E.mat');
if exist(fE,'file'), E = load(fE,'DE'); src = [src, num2cell(E.DE(1:2))]; end
% ids 13-14: yaw-constrained -E2 designs
fE2 = fullfile(simRoot,'..','results','yaw_ceiling','yaw_ceiling_designs.mat');
if exist(fE2,'file'), E2 = load(fE2,'DE2'); src = [src, num2cell(E2.DE2(1:2))]; end
% ids 15 PID-Budget-S, 16 A1-Commanded, 17-19 PID-E,
% AA-Budget-E, AA-HalfGain-E (estimator stack, iteration 0); ids 20-22 would hold
% a second -E iteration (designs_E_it1.mat); none exists, so ids 20-22 are the designs added below
RA = fullfile(simRoot,'..','results','analysis');
if exist(fullfile(RA,'designs_main.mat'),'file')
    Mm = load(fullfile(RA,'designs_main.mat')); Mm = Mm.M; P7 = load(fullfile(RA,'pidS.mat'),'PS'); P7 = P7.PS;
    src = [src, {keep(P7), keep(Mm.A1)}];
    E7 = load(fullfile(RA,'designs_E.mat'),'DE'); src = [src, arrayfun(@(x) keep(x), E7.DE, 'uni', 0)];
    if exist(fullfile(RA,'designs_E_it1.mat'),'file'), E8 = load(fullfile(RA,'designs_E_it1.mat'),'DE'); src = [src, arrayfun(@(x) keep(x), E8.DE, 'uni', 0)]; end
    % commanded-input designs with the ratio range
    % extended to 1.0, ids 20 A1-Commanded (ratio >= 1.0), 21 RO-HalfGain (ratio >= 1.0)
    % (no iteration-1 -E designs exist, so these follow id 19)
    if isfield(Mm, 'A1x') && isfield(Mm, 'ROx'), src = [src, {keep(Mm.A1x), keep(Mm.ROx)}]; end
    % id 22: PID-Budget with the b* set-point weights (setpoint_designs.mat, rule setpoint_rule.txt)
    if exist(fullfile(RA,'setpoint_designs.mat'),'file'), CB = load(fullfile(RA,'setpoint_designs.mat')); src = [src, {keep(CB.C(1).star)}]; end
end
names = cell(1,numel(src)); GT = zeros(37,numel(src)); KT = zeros(1,numel(src)); TH = zeros(2,numel(src));
for c = 1:numel(src)
    e = src{c}; names{c} = e.name; GT(:,c) = e.G(:); KT(c) = e.kind; TH(:,c) = [e.tauhat; e.ndhat];
end
save(fullfile(here,'sitl_gain_table.mat'),'GT','KT','TH','names');

mdl = 'sitl_ladrc_ctrl';
if bdIsLoaded(mdl), close_system(mdl,0); end
f = fullfile(here,[mdl '.slx']); if exist(f,'file'), delete(f); end
new_system(mdl); load_system(mdl);
set_param(mdl,'HardwareBoard','PX4 Host Target');
set_param(mdl,'SolverType','Fixed-step','Solver','FixedStepDiscrete','FixedStep','0.001','StopTime','inf');
mw = get_param(mdl,'ModelWorkspace');
assignin(mw,'GT',GT); assignin(mw,'KT',KT); assignin(mw,'TH',TH);

% --- uORB inputs
rd = {'lpos','VehicleLocalPosition','vehicle_local_position', 60;
      'att','VehicleAttitude','vehicle_attitude', 160;
      'avel','VehicleAngularVelocity','vehicle_angular_velocity', 260;
      'arm','ActuatorArmed','actuator_armed', 360;
      'odo','VehicleOdometry','vehicle_mocap_odometry', 600};
for i = 1:size(rd,1)
    b = [mdl '/' rd{i,1}];
    add_block('px4uORBlib/PX4 uORB Read', b, 'Position',[40 rd{i,4} 160 rd{i,4}+50]);
    set_param(b,'uORBMsg',rd{i,2},'uORBTopicName',rd{i,3},'SampleTime','0.001');
end
add_block('px4Utilslib/Read Parameter',[mdl '/pCtrl'],'Position',[40 440 160 480]);
set_param([mdl '/pCtrl'],'ParameterName','MC_ROLL_P','DataType','FLOAT','SampleTime','0.001');
add_block('px4Utilslib/Read Parameter',[mdl '/pScen'],'Position',[40 500 160 540]);
set_param([mdl '/pScen'],'ParameterName','MC_PITCH_P','DataType','FLOAT','SampleTime','0.001');
add_block('px4Utilslib/Read Parameter',[mdl '/pSrc'],'Position',[40 700 160 740]);
set_param([mdl '/pSrc'],'ParameterName','MC_YAWRATE_FF','DataType','FLOAT','SampleTime','0.001');

% --- bus selectors
sel = {'lpos','x,y,z,vx,vy,vz'; 'att','q'; 'avel','xyz'; 'arm','armed'; 'odo','position,velocity,q,angular_velocity'};
for i = 1:size(sel,1)
    bs = [mdl '/sel_' sel{i,1}];
    add_block('simulink/Signal Routing/Bus Selector', bs, 'Position',[230 rd{i,4}-10 240 rd{i,4}+70]);
    set_param(bs,'OutputSignals',sel{i,2});
    add_line(mdl,[sel{i,1} '/2'],['sel_' sel{i,1} '/1']);   % port 1 = Status, port 2 = Msg
    add_block('simulink/Sinks/Terminator',[mdl '/st_' sel{i,1}],'Position',[180 rd{i,4}-20 190 rd{i,4}-10]);
    add_line(mdl,[sel{i,1} '/1'],['st_' sel{i,1} '/1']);
end
add_block('simulink/Signal Routing/Mux',[mdl '/muxPos'],'Inputs','3','Position',[300 40 305 100]);
add_block('simulink/Signal Routing/Mux',[mdl '/muxVel'],'Inputs','3','Position',[300 110 305 170]);
for j = 1:3
    add_line(mdl,['sel_lpos/' num2str(j)],['muxPos/' num2str(j)]);
    add_line(mdl,['sel_lpos/' num2str(j+3)],['muxVel/' num2str(j)]);
end

% --- controller (MATLAB Function block)
fb = [mdl '/ctrl'];
add_block('simulink/User-Defined Functions/MATLAB Function', fb, 'Position',[420 150 600 750]);
rt = sfroot; ch = rt.find('-isa','Stateflow.EMChart','Path',fb);
ch.Script = sprintf(['function [u, dbg] = ctrl(pos, vel, q, gyro, armed, ctrlId, scen, opos, ovel, oq, orate, src)\n' ...
    '%%#codegen\n[u, dbg] = sitl_ctrl_step(pos, vel, q, gyro, armed, ctrlId, scen, GT, KT, TH, opos, ovel, oq, orate, src);\n']);
ch.SupportVariableSizing = 0;
for nm = {'GT','KT','TH'}
    d = Stateflow.Data(ch); d.Name = nm{1}; d.Scope = 'Parameter'; d.Tunable = false;
end
add_line(mdl,'muxPos/1','ctrl/1'); add_line(mdl,'muxVel/1','ctrl/2');
add_line(mdl,'sel_att/1','ctrl/3'); add_line(mdl,'sel_avel/1','ctrl/4');
add_line(mdl,'sel_arm/1','ctrl/5'); add_line(mdl,'pCtrl/1','ctrl/6'); add_line(mdl,'pScen/1','ctrl/7');
add_line(mdl,'sel_odo/1','ctrl/8'); add_line(mdl,'sel_odo/2','ctrl/9');
add_line(mdl,'sel_odo/3','ctrl/10'); add_line(mdl,'sel_odo/4','ctrl/11'); add_line(mdl,'pSrc/1','ctrl/12');
for nm = {'pCtrl','pScen','pSrc'}
    add_block('simulink/Sinks/Terminator',[mdl '/st_' nm{1}]); add_line(mdl,[nm{1} '/2'],['st_' nm{1} '/1']);
end

% --- actuator_motors output
add_block('px4uORBlib/PX4 uORB Message',[mdl '/motMsg'],'Position',[620 40 740 90]);
set_param([mdl '/motMsg'],'uORBMsg','ActuatorMotors','uORBTopicName','actuator_motors','SampleTime','0.001');
add_block('px4Utilslib/PX4 Timestamp',[mdl '/ts'],'Position',[620 110 700 140]);
set_param([mdl '/ts'],'SampleTime','0.001');
add_block('simulink/Signal Routing/Bus Assignment',[mdl '/busAsg'],'Position',[800 40 820 200]);
set_param([mdl '/busAsg'],'AssignedSignals','timestamp,timestamp_sample,control');
add_line(mdl,'motMsg/1','busAsg/1'); add_line(mdl,'ts/1','busAsg/2'); add_line(mdl,'ts/1','busAsg/3');
add_line(mdl,'ctrl/1','busAsg/4');
add_block('px4uORBlib/PX4 uORB Write',[mdl '/motWrite'],'Position',[880 90 1000 150]);
set_param([mdl '/motWrite'],'uORBMsg','ActuatorMotors','uORBTopicName','actuator_motors');
add_line(mdl,'busAsg/1','motWrite/1');
% --- debug_array (measurement and virtual controls at every step, logging only)
add_block('px4uORBlib/PX4 uORB Message',[mdl '/dbgMsg'],'Position',[620 300 740 350]);
set_param([mdl '/dbgMsg'],'uORBMsg','DebugArray','uORBTopicName','debug_array','SampleTime','0.001');
add_block('simulink/Signal Routing/Bus Assignment',[mdl '/dbgAsg'],'Position',[800 300 820 420]);
set_param([mdl '/dbgAsg'],'AssignedSignals','timestamp,data');
add_line(mdl,'dbgMsg/1','dbgAsg/1'); add_line(mdl,'ts/1','dbgAsg/2'); add_line(mdl,'ctrl/2','dbgAsg/3');
add_block('px4uORBlib/PX4 uORB Write',[mdl '/dbgWrite'],'Position',[880 330 1000 390]);
set_param([mdl '/dbgWrite'],'uORBMsg','DebugArray','uORBTopicName','debug_array');
add_line(mdl,'dbgAsg/1','dbgWrite/1');
save_system(mdl, f);
end

function s = keep(c)
s = struct('name', c.name, 'kind', c.kind, 'G', c.G, 'tauhat', c.tauhat, 'ndhat', c.ndhat);
end
