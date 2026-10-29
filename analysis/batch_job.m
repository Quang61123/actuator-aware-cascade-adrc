function batch_job(name, nw)
if nargin < 2, nw = 11; end
%BATCH_JOB Batch runner for long computations: evaluates the command
%   in <work dir>\<name>.cmd.txt and writes
%   <name>.done (JOB_OK or JOB_FAIL with the error report).
W = px4sitl_work_dir(); here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here,'..','MATLAB_simulation'), fullfile(here,'..','estimator_noise'), fullfile(here,'..','yaw_ceiling'));
done = fullfile(W, ['' name '.done']);
try
    try, if isempty(gcp('nocreate')), cl = parcluster('Processes'); cl.NumWorkers = nw; parpool(cl, nw); end, catch ep, fprintf('NO_POOL (%s): parfor runs serially\n', ep.message); end
    eval(fileread(fullfile(W, ['' name '.cmd.txt'])));
    f = fopen(done,'w'); fprintf(f,'JOB_OK\n'); fclose(f);
catch e
    f = fopen(done,'w'); fprintf(f,'JOB_FAIL %s\n', getReport(e,'extended','hyperlinks','off')); fclose(f);
end
end
