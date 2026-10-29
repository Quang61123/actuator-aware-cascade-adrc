function audit_design()
P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));allC=[D.C D.Strict];T=table();GT=table();F=table();chs=[1 3 4 6];for j=1:numel(allC)
 c=allC(j);coCap=D.cap;if j>7,coCap=D.strict;end;zz=cascade_stats(c,D.S);assert(zz(4)<=coCap*(1+1e-9));F=[F;table(string(c.name),coCap,zz(4),zz(1),'VariableNames',{'Ctrl','CascadeCap','PredSigmaDU','SigmaE'})];
 for h=1:4,g=c.G(6*(chs(h)-1)+(1:6))';gi=[];if h==1,gi=c.G(19:24)';end;z=channel_stats(c.kind,g,h,D.S,gi,[c.tauhat c.ndhat],D.S.act,true,P.nfreq);lim=D.dub(h);if j==1||j==7,lim=Inf;end;fprintf('AUDIT %s class=%d z=%s limDU=%.9g\n',c.name,h,mat2str(z,9),lim);assert(all(isfinite(z))&&z(4)<=2&&z(5)<=10&&z(2)<=D.S.ub(h)*(1+1e-8)&&z(3)<=lim*(1+1e-8));T=[T;table(string(c.name),h,z(1),z(2),z(3),z(4),z(5),z(6),'VariableNames',{'Ctrl','Class','SigmaE','SigmaU','SigmaDU','Ms','Overshoot','Rho'})];GT=[GT;table(string(c.name),h,g(1),g(2),g(3),g(4),g(5),g(6),'VariableNames',{'Ctrl','Class','G1','G2','G3','G4','G5','G6'})];end
end
writetable(T,fullfile(P.results,'design_audit.csv'));writetable(GT,fullfile(P.results,'gains.csv'));writetable(F,fullfile(P.results,'cascade_design_audit.csv'));
end
