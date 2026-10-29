function guard=simulink_cache()
% Keep transient compilation files on the local disk; restore session settings.
cfg=Simulink.fileGenControl('getConfig');previousCache=cfg.CacheFolder;previousCode=cfg.CodeGenFolder;
local=fullfile(tempdir,'simulink_cache');
Simulink.fileGenControl('set','CacheFolder',fullfile(local,'cache'),...
 'CodeGenFolder',fullfile(local,'codegen'),'createDir',true);
guard=onCleanup(@()Simulink.fileGenControl('set','CacheFolder',previousCache,'CodeGenFolder',previousCode));
end
