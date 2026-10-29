function P=sim_config()
%SIM_CONFIG Simulation configuration (paths, seeds, run counts, grids).
% Paths are resolved relative to this file: results in ../results,
% figures in ../figures.
P.root=fileparts(mfilename('fullpath'));
P.results=fullfile(P.root,'..','results');
P.figures=fullfile(P.root,'..','figures');
P.workers=2;P.Ts=.001;P.duration=40;P.window=[10 40];
P.nominalSeed=91000;P.stressSeed=93000;P.gridSeed=92000;
P.mcRuns=30;P.stressRuns=10;P.gridRuns=10;
P.taus=[.015 .03 .045 .06 .075 .09 .105 .12];P.delays=[0 5 10];
P.nfreq=4000;P.oscThreshold=20;
end
