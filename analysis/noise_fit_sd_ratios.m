function T = noise_fit_sd_ratios()
%NOISE_FIT_SD_RATIOS Converts the held-out increment-variance errors of the 60-s per-design
%   noise models (estimator_noise/est_noise_validation.csv, rows with a design id) into standard-deviation
%   ratios: SD error = sqrt(1+e/100)-1. Rate outputs (e_*dot) and the other outputs reported separately.
%   Writes noise_fit_sd_ratios.csv.
here = fileparts(mfilename('fullpath')); RP = fullfile(here,'..','results');
V = readtable(fullfile(RP,'estimator_noise','est_noise_validation.csv'),'TextType','string'); V = V(~isnan(V.Designs),:);
V.SDErrPct = 100*(sqrt(1+V.ValidIncVarErrPct/100)-1); rate = endsWith(V.Output,"dot");
oth = V(~rate,:); [~,k] = max(abs(oth.ValidIncVarErrPct)); o2 = oth; o2(k,:) = [];
T = table(["rates";"others (excl. max)";"max other"], [max(abs(V.ValidIncVarErrPct(rate)));max(abs(o2.ValidIncVarErrPct));abs(oth.ValidIncVarErrPct(k))], ...
  [max(abs(V.SDErrPct(rate)));max(abs(o2.SDErrPct));oth.SDErrPct(k)], [min(V.SDErrPct(rate));min(o2.SDErrPct);NaN], [max(V.SDErrPct(rate));max(o2.SDErrPct);NaN], ...
  'VariableNames',{'Group','MaxAbsVarErrPct','MaxAbsSDErrPct','MinSDErrPct','MaxSDErrPct'});
fprintf('max other: %s %s design %d\n', oth.Class(k), oth.Output(k), oth.Designs(k));
writetable(T, fullfile(RP,'analysis','noise_fit_sd_ratios.csv')); disp(T)
end