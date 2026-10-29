function linear_checks()
P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));C=D.C;S=D.S;
% E2: controlled roll-channel bandwidth sweeps (no retuning during each curve).
r=logspace(log10(1.5),log10(60),50);Sweep=table();
for wc=[3 6 10],for kind=[4 7],for rr=r
 g=gains_of(kind,rr,wc);[ms,rho]=chan_ms(kind,g,S.Ts,[],S.act,[],S.act,P.nfreq);Sweep=[Sweep;table(wc,kind,rr,ms,rho,'VariableNames',{'WC','Kind','Ratio','Ms','Rho'})];
end,end,end
writetable(Sweep,fullfile(P.results,'sensitivity_ratio.csv'));
% E6: full channel-class limits for fixed nominal gains, never robust retuned.
Mismatch=table();ch=[1 3 4 6];
for c=[1 2 4],for h=1:4
 g=C(c).G(6*(ch(h)-1)+(1:6))';gi=[];if h==1,gi=C(c).G(19:24)';end;oa=[C(c).tauhat C(c).ndhat];
 for tau=.01:.0025:.15
  [ms,rho]=chan_ms(C(c).kind,g,S.Ts,gi,[tau 5],[],oa,P.nfreq);Mismatch=[Mismatch;table(string(C(c).name),h,string('lag'),tau*1000,ms,rho,'VariableNames',{'Ctrl','Class','Axis','Value','Ms','Rho'})];
 end
 for nd=0:30
  [ms,rho]=chan_ms(C(c).kind,g,S.Ts,gi,[.03 nd],[],oa,P.nfreq);Mismatch=[Mismatch;table(string(C(c).name),h,string('delay'),nd,ms,rho,'VariableNames',{'Ctrl','Class','Axis','Value','Ms','Rho'})];
 end
end,end
writetable(Mismatch,fullfile(P.results,'sensitivity_mismatch.csv'));
% Freshly audit every inherited curve gain under exactly its stated local ceiling.
A=table();for j=1:numel(D.Curve),cc=D.Curve(j).Controller;for h=1:4
 g=cc.G(6*(ch(h)-1)+(1:6))';gi=[];if h==1,gi=cc.G(19:24)';end;z=channel_stats(cc.kind,g,h,S,gi,[cc.tauhat cc.ndhat],S.act,true,P.nfreq);assert(all(isfinite(z))&&z(4)<=2&&z(5)<=10&&z(2)<=S.ub(h)*(1+1e-8)&&z(3)<=D.Curve(j).Gamma*D.dub(h)*(1+1e-8));
 A=[A;table(string(D.Curve(j).Family),D.Curve(j).Gamma,h,z(1),z(2),z(3),z(4),z(5),'VariableNames',{'Family','Gamma','Class','SigmaE','SigmaU','SigmaDU','Ms','Overshoot'})];
end,end
writetable(A,fullfile(P.results,'budget_design_audit.csv'));save(fullfile(P.results,'linear_checks.mat'),'Sweep','Mismatch','A','C','S','P');fprintf('LINEAR_CHECKS_PASS\n');
end
