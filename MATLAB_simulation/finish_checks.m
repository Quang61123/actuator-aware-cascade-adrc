function finish_checks()
here=fileparts(mfilename('fullpath'));cd(here);addpath(here,fullfile(here,'..','Simulink'));validate_simulink_model();validate_mex40();validate_outputs();validate_cascade();make_figures();export_results();fprintf('CHECKS_PASS\n');
end
