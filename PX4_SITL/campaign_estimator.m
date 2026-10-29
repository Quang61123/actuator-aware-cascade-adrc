%CAMPAIGN_ESTIMATOR Estimator-stack SITL runs (run_sitl_estimator); parallel: set SITL_WORKER / SITL_NWORKERS / SITL_PARTS.
here = fileparts(mfilename('fullpath')); addpath(here);
w = str2double(getenv('SITL_WORKER')); nw = str2double(getenv('SITL_NWORKERS'));
if isnan(w), w = 0; nw = 1; end
parts = getenv('SITL_PARTS'); if isempty(parts), parts = 'HS'; end
lg = fullfile(px4sitl_work_dir(), sprintf('est_w%d.log', w));
try
    run_sitl_estimator(parts, w, nw);
    fid = fopen(lg,'a'); fprintf(fid,'EST_DONE worker %d\n', w); fclose(fid);
catch e
    fid = fopen(lg,'a'); fprintf(fid,'EST_FAIL %s\n', getReport(e,'extended')); fclose(fid);
end
