function [W, Wwsl] = px4sitl_work_dir()
%PX4SITL_WORK_DIR Local working folder of the PX4 SITL setup (logs, ulogs, build, launch scripts).
%   Set the environment variable PX4SITL_WORK to override; default %USERPROFILE%\px4sitl_work.
%   Wwsl is the same folder as seen from WSL (/mnt/<drive>/...).
W = getenv('PX4SITL_WORK'); if isempty(W), W = fullfile(getenv('USERPROFILE'), 'px4sitl_work'); end
Wwsl = ['/mnt/' lower(W(1)) strrep(W(3:end), '\', '/')];
end