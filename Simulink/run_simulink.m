function out=run_simulink(scenario,controller,duration)
% Run the model using shared MATLAB simulation functions.
if nargin<1,scenario=6;end
if nargin<2,controller='AA-Budget';end
here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'..','MATLAB_simulation'));
cacheGuard=simulink_cache(); %#ok<NASGU>
model=fullfile(here,'control_system.slx');
if ~isfile(model),create_simulink_model;end
if nargin<3,setup_simulink(scenario,controller);else,setup_simulink(scenario,controller,duration);end
load_system(model);out=sim('control_system');
end
