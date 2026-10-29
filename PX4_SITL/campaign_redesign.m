%CAMPAIGN_REDESIGN Redesign SITL runs from run_plan.csv (run_sitl_redesign); set SITL_WORKER / SITL_NWORKERS / SITL_PHASE / SITL_EF.
here = fileparts(mfilename('fullpath')); addpath(here);
w = str2double(getenv('SITL_WORKER')); nw = str2double(getenv('SITL_NWORKERS')); if isnan(w), w = 0; nw = 1; end
ph = getenv('SITL_PHASE'); ef = str2num(getenv('SITL_EF')); if isempty(ef), ef = [17 18 19]; end %#ok<ST2NM>
lg = fullfile(px4sitl_work_dir(), sprintf('w%d.log', w));
try
    wt = getenv('SITL_WAITFOR'); % wait until this worker finished the given phase
    while ~isempty(wt) && ~(exist(lg,'file') && contains(fileread(lg), sprintf('RUN_DONE phase %s worker %d', wt, w))), pause(20); end
    run_sitl_redesign(ph, w, nw, ef);
catch e
    fid = fopen(lg,'a'); fprintf(fid,'RUN_FAIL %s\n', getReport(e,'extended')); fclose(fid);
end
