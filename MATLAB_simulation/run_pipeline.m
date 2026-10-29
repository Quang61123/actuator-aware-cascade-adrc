function run_pipeline()
here=fileparts(mfilename('fullpath'));cd(here);addpath(here);run_design();run_compute();finish_checks();fprintf('PIPELINE_COMPLETE\n');
end
