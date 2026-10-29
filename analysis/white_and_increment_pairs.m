function T = white_and_increment_pairs()
%WHITE_AND_INCREMENT_PAIRS (1) Paired rotor-command increment J of AA-Budget-E / AA-HalfGain-E versus PID-E in SITL
%   (S6 EFS 30 pairs, S5 EF5 10 pairs, from sitl_all_runs.csv); (2) white-noise predicted roll/pitch
%   increment of the MATLAB designs over the measured 300-s held-out SITL value (validation_3level.csv).
%   Writes J_vsPIDE.csv and white_over_measured.csv.
here = fileparts(mfilename('fullpath')); RA = fullfile(here,'..','results','analysis');
A = readtable(fullfile(RA,'sitl_all_runs.csv'),'TextType','string'); T = table();
for g = ["EFS","EF5"], P = A(A.Group==g & A.Ctrl==17,:);
  for c = [18 19], Q = A(A.Group==g & A.Ctrl==c,:); [~,ia,ib] = intersect(P.Seed,Q.Seed); d = Q.J_unq(ib)-P.J_unq(ia); n = numel(d);
    T = [T; table(g, c, n, mean(P.J_unq(ia)), mean(Q.J_unq(ib)), mean(d), tinv(0.975,n-1)*std(d)/sqrt(n), sum(d<0), 'VariableNames',{'Group','Ctrl','N','J_PIDE','J','Diff','CI','PairsLower'})]; end, end
writetable(T, fullfile(RA,'J_vsPIDE.csv')); disp(T)
V = readtable(fullfile(RA,'validation_3level.csv'),'TextType','string');
w = {"PID-Budget",1.380; "AA-Budget",0.937; "AA-HalfGain",1.023; "PID-Budget-S",0.5917; "AA-Budget-Strict",0.496; "AA-HalfGain-Strict",0.521}; W = table();
for k = 1:size(w,1), m = max(V.Measured_heldout(V.Design==w{k,1} & ismember(V.Quantity,["sdu_phi","sdu_theta"])));
  W = [W; table(w{k,1}, w{k,2}, m, w{k,2}/m, 'VariableNames',{'Design','WhitePred','Measured300s','Factor'})]; end
writetable(W, fullfile(RA,'white_over_measured.csv')); disp(W)
end