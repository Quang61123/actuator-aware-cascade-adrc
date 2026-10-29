function run_compute()
here=fileparts(mfilename('fullpath'));cd(here);addpath(here);P=sim_config();if license('test','MATLAB_Coder'),build_extra_mex();else,warning('sim:NoCoder','Coder unavailable: MATLAB fallback will be slower.');end;if isempty(gcp('nocreate')),parpool('Processes',P.workers);end;mc_trials();linear_checks();p=gcp('nocreate');if ~isempty(p),delete(p);end;numerical_checks();fprintf('COMPUTE_COMPLETE\n');
end
