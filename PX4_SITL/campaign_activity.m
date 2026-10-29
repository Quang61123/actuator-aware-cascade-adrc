%CAMPAIGN_ACTIVITY Activity-check SITL runs (run_sitl_activity); parallel: set SITL_WORKER / SITL_NWORKERS.
here = fileparts(mfilename('fullpath')); addpath(here);
w = str2double(getenv('SITL_WORKER')); nw = str2double(getenv('SITL_NWORKERS'));
if isnan(w), w = 0; nw = 1; end
parts = getenv('SITL_PARTS'); if isempty(parts), parts = 'SH'; end
lg = fullfile(px4sitl_work_dir(), sprintf('act_w%d.log', w));
try
    run_sitl_activity(parts, w, nw);
    fid = fopen(lg,'a'); fprintf(fid,'ACT_DONE worker %d\n', w); fclose(fid);
catch e
    fid = fopen(lg,'a'); fprintf(fid,'ACT_FAIL %s\n', getReport(e,'extended')); fclose(fid);
end
