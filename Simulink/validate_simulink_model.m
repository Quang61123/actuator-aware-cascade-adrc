function V=validate_simulink_model()
here=fileparts(mfilename('fullpath'));root=fullfile(here,'..','MATLAB_simulation');addpath(here,root);cacheGuard=simulink_cache();P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));C=D.C;S=D.S;create_simulink_model;model='control_system';V=table();
for c=1:4
 pv=quad_params();pc=quad_params([],C(c).tauhat,C(c).ndhat);E=make_exo(6,P.nominalSeed+1,40,S.Ts,S);N=40000;L=run_quad_sim(C(c).kind,C(c).G,6,pv,pc,40,S.Ts,E,S.sig);
 in=Simulink.SimulationInput(model);in=in.setVariable('Ts',S.Ts).setVariable('Tend',40-S.Ts).setVariable('scen',6).setVariable('kind',C(c).kind).setVariable('G',C(c).G).setVariable('pv',pc).setVariable('pvt',pv).setVariable('sig',S.sig).setVariable('exo_ts',[(0:N-1)'*S.Ts,[E.noise;E.bias;E.wind]']);out=sim(in);x=squeeze(out.sl_x);if size(x,1)==16,x=x';end;assert(isequal(size(x),size(L.x))&&all(isfinite(x),'all'));delta=abs(x-L.x);assert(max(delta,[],'all')<1e-8);V=[V;table(string(C(c).name),N,max(delta(:,1:3),[],'all'),max(delta(:,7:9),[],'all'),max(delta(:,13:16),[],'all'),'VariableNames',{'Ctrl','Samples','MaxPositionDiff','MaxAttitudeDiff','MaxRotorDiff'})];fprintf('SIMULINK %s max=%.3g\n',C(c).name,max(delta,[],'all'));
end
writetable(V,fullfile(P.results,'simulink.csv'));save(fullfile(P.results,'simulink.mat'),'V');fprintf('SIMULINK_PASS\n');
end
