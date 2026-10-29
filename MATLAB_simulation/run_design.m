function run_design()
here=fileparts(mfilename('fullpath'));cd(here);addpath(here);P=sim_config();if isempty(gcp('nocreate')),parpool('Processes',P.workers);end;primary_design();p=gcp('nocreate');if ~isempty(p),delete(p);end
end
