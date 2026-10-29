%DO_REPORT_MAIN Run sitl_report_main in a separate MATLAB process.
here = fileparts(mfilename('fullpath')); addpath(here);
lg = fullfile(px4sitl_work_dir(),'report_main.log');
try
    Rep = sitl_report_main(); %#ok<NASGU>
    fid = fopen(lg,'a'); fprintf(fid,'REPORT_OK\n'); fclose(fid);
catch e
    fid = fopen(lg,'a'); fprintf(fid,'REPORT_FAIL %s\n', getReport(e,'extended')); fclose(fid);
end
