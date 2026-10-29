function build_mex()
here=fileparts(mfilename('fullpath'));cd(here);P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));S=D.S;pv=quad_params();E=make_exo(6,91001,40,S.Ts,S);nd=5*ones(40000,1);
assert(license('test','MATLAB_Coder'),'MATLAB Coder unavailable');disp(mex.getCompilerConfigurations('C++','Selected'));cfg=coder.config('mex');cfg.GenerateReport=false;
codegen('run_quad_stress','-config',cfg,'-args',{D.C(2).kind,D.C(2).G,6,pv,pv,40,S.Ts,E,S.sig,nd},'-o','run_quad_stress_mex');
A=table();for c=1:4
 pc=quad_params([],D.C(c).tauhat,D.C(c).ndhat);a=run_quad_stress(D.C(c).kind,D.C(c).G,6,pv,pc,40,S.Ts,E,S.sig,nd);b=run_quad_stress_mex(D.C(c).kind,D.C(c).G,6,pv,pc,40,S.Ts,E,S.sig,nd);delta=max(abs(a.x-b.x),[],'all');assert(delta<1e-6);A=[A;table(string(D.C(c).name),delta,'VariableNames',{'Ctrl','MaxStateDifference'})];
end
writetable(A,fullfile(P.results,'mex_equivalence.csv'));save(fullfile(P.results,'mex_equivalence.mat'),'A');fprintf('MEX_EQUIVALENCE_PASS\n');
end
