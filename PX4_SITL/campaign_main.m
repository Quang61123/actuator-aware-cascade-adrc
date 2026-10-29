%CAMPAIGN_MAIN Run the SITL campaign. Parallel use: set the
%   environment variables SITL_WORKER (0..n-1) and SITL_NWORKERS before
%   matlab -batch "run('campaign_main.m')"; worker w uses PX4 instance w.
here = fileparts(mfilename('fullpath')); addpath(here);
w = str2double(getenv('SITL_WORKER')); nw = str2double(getenv('SITL_NWORKERS'));
if isnan(w), w = 0; nw = 1; end
parts = 'ABCDEF';
lg = fullfile(px4sitl_work_dir(), sprintf('campaign_main_w%d.log', w));
try
    run_sitl_campaign(parts, w, nw);
    fid = fopen(lg,'a'); fprintf(fid,'CAMPAIGN_MAIN_DONE worker %d\n', w); fclose(fid);
catch e
    fid = fopen(lg,'a'); fprintf(fid,'CAMPAIGN_MAIN_FAIL %s\n', getReport(e,'extended')); fclose(fid);
end
