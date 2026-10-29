function mc_paired()
P=sim_config();D=load(fullfile(P.results,'trials.mat'));T=D.T;A=table();pairs={'AA-Budget','PID-Budget';'AA-HalfGain','AA-Budget';'A1-Commanded','AA-Budget';'A2-Unaligned','AA-Budget';'A3-NoIncrement','AA-Budget'};
for j=1:size(pairs,1)
 a=T(ismember(string(T.Experiment),["E1","E3"])&string(T.Ctrl)==string(pairs{j,1}),:);b=T(string(T.Experiment)=="E1"&string(T.Ctrl)==string(pairs{j,2}),:);[a,ia]=sortrows(a,'Seed');[b,ib]=sortrows(b,'Seed');assert(isequal(a.Seed,b.Seed));ok=~a.Failed&~b.Failed;v=a.RMSE_mm(ok)-b.RMSE_mm(ok);u=a.Increment(ok)-b.Increment(ok);n=numel(v);ci=tinv(.975,n-1);A=[A;table(string(pairs{j,1}),string(pairs{j,2}),n,mean(v),ci*std(v)/sqrt(n),mean(u),ci*std(u)/sqrt(n),'VariableNames',{'A','B','N','RMSEDiff_mm','RMSEDiffCI95_mm','IncrementDiff','IncrementDiffCI95'})];
end
writetable(A,fullfile(P.results,'paired.csv'));save(fullfile(P.results,'paired.mat'),'A');disp(A);fprintf('PAIRED_PASS\n');
end
