function mc_summarize()
P=sim_config();D=load(fullfile(P.results,'trials.mat'));T=D.T;gammaToken=repmat("none",height(T),1);sel=isfinite(T.Gamma);gammaToken(sel)=string(T.Gamma(sel));key=string(T.Experiment)+"|"+string(T.Ctrl)+"|"+string(T.Tau)+"|"+string(T.Delay)+"|"+string(T.Noise)+"|"+gammaToken;groups=unique(key,'stable');U=table();
for j=1:numel(groups)
 a=T(key==groups(j),:);ok=~a.Failed;v=a.RMSE_mm(ok);n=numel(v);ci=NaN;if n>1,ci=tinv(.975,n-1)*std(v)/sqrt(n);end
 U=[U;table(string(a.Experiment(1)),string(a.Ctrl(1)),a.Tau(1)*1000,a.Delay(1),a.Noise(1),a.Gamma(1),height(a),sum(a.Failed),sum(a.FluctuationThreshold),mean(v),std(v),ci,mean(a.Increment(ok)),mean(a.RotorStd(ok)),mean(a.RotorClip_pct(ok)),mean(a.VirtualClip_pct(ok)),'VariableNames',{'Experiment','Ctrl','Tau_ms','Delay_ms','NoiseScale','Gamma','Runs','Failures','ThresholdExceedances','RMSE_mean_mm','RMSE_std_mm','RMSE_ci95_mm','Increment','RotorStd','RotorClip_pct','VirtualClip_pct'})];
end
U.MaxError_mean_mm=zeros(height(U),1);U.MaxError_std_mm=U.MaxError_mean_mm;U.MaxError_ci95_mm=U.MaxError_mean_mm;for j=1:numel(groups),a=T(key==groups(j)&~T.Failed,:);v=a.MaxError_mm;n=numel(v);U.MaxError_mean_mm(j)=mean(v);U.MaxError_std_mm(j)=std(v);if n>1,U.MaxError_ci95_mm(j)=tinv(.975,n-1)*std(v)/sqrt(n);else,U.MaxError_ci95_mm(j)=NaN;end,end
writetable(U,fullfile(P.results,'summary.csv'));save(fullfile(P.results,'summary.mat'),'U');
end
