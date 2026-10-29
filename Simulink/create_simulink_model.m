here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'..','MATLAB_simulation'));
%% CREATE_SIMULINK_MODEL  Builds control_system.slx programmatically.
%  Cascade controller and quadrotor plant with rotor dynamics (16-state plant), delayed sensor buffer,
%  exogenous noise / gyro bias / wind signals from the workspace.
% Initialize numeric plant parameters and all variables needed by the model.
setup_simulink;
mdl = 'control_system';
if bdIsLoaded(mdl), close_system(mdl, 0); end
% save_system overwrites only the destination model in this folder.
new_system(mdl);
set_param(mdl, 'Solver','ode4', 'SolverType','Fixed-step', 'FixedStep','Ts', 'StopTime','Tend', ...
    'SaveOutput','off','SaveTime','off','SignalLogging','off');
addb = @(lib, name, pos) add_block(lib, [mdl '/' name], 'Position', pos);
fcn = @(name, pos, code) add_fcn(mdl, name, pos, code);
konst = @(name, val, pos) set_param(add_block('simulink/Sources/Constant', [mdl '/' name], 'Position', pos), 'Value', val);
% sources
addb('simulink/Sources/Clock', 'Clock', [20 30 50 60]);
konst('scen', 'scen', [20 90 60 110]);
fcn('Reference', [100 20 200 80], sprintf('function r = fcn(t, scen)\n%%#codegen\nr = reference_gen(t, scen);\n'));
fcn('Disturbance', [560 -120 680 -40], sprintf('function [Fd, Td] = fcn(t, scen, exo)\n%%#codegen\n[Fd, Td] = disturbance_gen(t, scen);\nFd = Fd + exo(16:18);\n'));
addb('simulink/Sources/From Workspace', 'Exogenous (noise, bias, wind)', [930 330 1030 360]);
set_param([mdl '/Exogenous (noise, bias, wind)'], 'VariableName','exo_ts','SampleTime','Ts','Interpolate','off', ...
    'OutputAfterFinalValue','Holding final value','Orientation','left');
% controller
fcn('Cascade Controller', [300 20 440 200], sprintf(['function [T, tau, att_ref, csn] = fcn(meas, r, cs, kind, G, pv, Ts)\n%%#codegen\n' ...
    '[T, tau, att_ref, csn] = ctrl_law(kind, G, cs, meas, r, pv, Ts);\n']));
addb('simulink/Discrete/Unit Delay', 'Controller State', [345 240 395 280]);
set_param([mdl '/Controller State'], 'InitialCondition','zeros(96,1)', 'SampleTime','Ts', 'Orientation','left');
konst('kind','kind',[230 100 260 115]); konst('Gains','G',[230 125 260 140]);
konst('pv nominal','pv',[230 150 260 165]); konst('Ts','Ts',[230 175 260 190]);
% mixer + plant
fcn('Rotor Mixer + Saturation', [500 40 610 110], sprintf('function wcmd = fcn(T, tau, pvt)\n%%#codegen\n[wcmd, ~, ~] = quad_mixer(T, tau, pvt);\n'));
konst('pv true','pvt',[430 -20 470 0]);
fcn('Quadrotor + Rotor Dynamics', [720 30 860 150], sprintf('function xdot = fcn(x, wcmd, Fd, Td, pvt)\n%%#codegen\nxdot = quad_dynamics(x, wcmd, Fd, Td, pvt);\n'));
addb('simulink/Continuous/Integrator', 'States x', [900 70 935 110]);
set_param([mdl '/States x'], 'InitialCondition','[zeros(12,1); sqrt(pvt(1)*pvt(2)/(4*pvt(4)))*ones(4,1)]');
% sensors with delay buffer
fcn('Sensors + Delay', [700 270 840 380], sprintf(['function [meas, bufn, initn] = fcn(x, exo, sig, buf, init, pvt)\n%%#codegen\n' ...
    '[meas, bufn, initn] = sensor_buffer(x, exo(1:12), sig, exo(13:15), buf, init, pvt(17));\n']));
set_param([mdl '/Sensors + Delay'], 'Orientation','left');
addb('simulink/Discrete/Unit Delay', 'Sensor Buffer', [760 410 810 440]); set_param([mdl '/Sensor Buffer'],'InitialCondition','zeros(12,11)','SampleTime','Ts');
addb('simulink/Discrete/Unit Delay', 'Init Flag', [760 460 810 490]); set_param([mdl '/Init Flag'],'InitialCondition','0','SampleTime','Ts');
konst('sigma','sig',[950 390 990 410]); konst('pv true 2','pvt',[950 430 990 450]);
% logging
addb('simulink/Sinks/Scope', 'State Scope', [1010 60 1040 100]);
addb('simulink/Sinks/To Workspace', 'log x', [1000 140 1060 160]);
set_param([mdl '/log x'],'VariableName','sl_x','SaveFormat','Array','SampleTime','Ts','MaxDataPoints','inf');
addb('simulink/Sinks/To Workspace', 'log wcmd', [700 -30 760 -10]);
set_param([mdl '/log wcmd'],'VariableName','sl_w','SaveFormat','Array','SampleTime','Ts','MaxDataPoints','inf');
for b = {'Reference','Disturbance','Cascade Controller','Rotor Mixer + Saturation','Sensors + Delay'}
    set_param([mdl '/' b{1}], 'SystemSampleTime', 'Ts');
end
L = @(a, b) add_line(mdl, a, b, 'autorouting', 'smart');
L('Clock/1','Reference/1'); L('scen/1','Reference/2');
L('Clock/1','Disturbance/1'); L('scen/1','Disturbance/2'); L('Exogenous (noise, bias, wind)/1','Disturbance/3');
L('Sensors + Delay/1','Cascade Controller/1'); L('Reference/1','Cascade Controller/2');
L('Controller State/1','Cascade Controller/3'); L('kind/1','Cascade Controller/4'); L('Gains/1','Cascade Controller/5');
L('pv nominal/1','Cascade Controller/6'); L('Ts/1','Cascade Controller/7'); L('Cascade Controller/4','Controller State/1');
L('Cascade Controller/1','Rotor Mixer + Saturation/1'); L('Cascade Controller/2','Rotor Mixer + Saturation/2'); L('pv true/1','Rotor Mixer + Saturation/3');
L('Rotor Mixer + Saturation/1','Quadrotor + Rotor Dynamics/2'); L('Rotor Mixer + Saturation/1','log wcmd/1');
L('Disturbance/1','Quadrotor + Rotor Dynamics/3'); L('Disturbance/2','Quadrotor + Rotor Dynamics/4'); L('pv true/1','Quadrotor + Rotor Dynamics/5');
L('Quadrotor + Rotor Dynamics/1','States x/1'); L('States x/1','Quadrotor + Rotor Dynamics/1');
L('States x/1','State Scope/1'); L('States x/1','log x/1');
L('States x/1','Sensors + Delay/1'); L('Exogenous (noise, bias, wind)/1','Sensors + Delay/2'); L('sigma/1','Sensors + Delay/3');
L('Sensor Buffer/1','Sensors + Delay/4'); L('Init Flag/1','Sensors + Delay/5'); L('pv true 2/1','Sensors + Delay/6');
L('Sensors + Delay/2','Sensor Buffer/1'); L('Sensors + Delay/3','Init Flag/1');
save_system(mdl,fullfile(here,[mdl '.slx'])); fprintf('Created %s.slx\n', mdl);

function add_fcn(mdl, name, pos, code)
add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/' name], 'Position', pos);
rt = sfroot; ch = rt.find('-isa','Stateflow.EMChart','Path',[mdl '/' name]); ch.Script = code;
end



