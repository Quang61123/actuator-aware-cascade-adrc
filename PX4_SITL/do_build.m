%DO_BUILD Generate code for the controller model and build PX4 SITL
%   with it (matlab -batch "run('do_build.m')"). This replaces the
%   px4_simulink_app module of the shared PX4 build in WSL.
W = fullfile(px4sitl_work_dir(),'build');
if ~exist(W,'dir'), mkdir(W); end
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'));
cd(W);
t0 = tic;
try
    mdl = create_sitl_controller();
    slbuild(mdl);
    fid = fopen(fullfile(px4sitl_work_dir(),'build.log'),'a'); fprintf(fid,'BUILD_OK %.1f s\n', toc(t0)); fclose(fid);
catch e
    fid = fopen(fullfile(px4sitl_work_dir(),'build.log'),'a'); fprintf(fid,'BUILD_FAIL %s\n', getReport(e,'extended')); fclose(fid);
end
