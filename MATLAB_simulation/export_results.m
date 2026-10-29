function export_results()
%EXPORT_RESULTS Adds the MeasuredPass column (measured cascade increment within the cap) to cascade_noise.csv.
P=sim_config();
K=readtable(fullfile(P.results,'cascade_noise.csv'));K.MeasuredPass=max(K.MeasuredSigmaDU,K.PitchMeasuredSigmaDU)<=K.CascadeCap*(1+1e-9);writetable(K,fullfile(P.results,'cascade_noise.csv'));fprintf('RESULTS_EXPORT_PASS\n');
end
