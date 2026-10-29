function build_extra_mex()
P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));S=D.S;pv=quad_params();cfg=coder.config('mex');cfg.GenerateReport=false;E=make_exo(5,91001,45,S.Ts,S);nd=5*ones(45000,1);codegen('run_quad_stress','-config',cfg,'-args',{D.C(2).kind,D.C(2).G,5,pv,pv,45,S.Ts,E,S.sig,nd},'-o','run_quad_s5_mex');
a=run_quad_stress(D.C(2).kind,D.C(2).G,5,pv,pv,45,S.Ts,E,S.sig,nd);b=run_quad_s5_mex(D.C(2).kind,D.C(2).G,5,pv,pv,45,S.Ts,E,S.sig,nd);assert(isequaln(a.x,b.x)&&isequaln(a.wcmd,b.wcmd));Q=table("S5",max(abs(a.x-b.x),[],'all'),'VariableNames',{'Case','MaxStateDifference'});
E=make_exo(2,96001,60,S.Ts,S);E.bias(:)=0;E.wind(:)=0;nd=5*ones(60000,1);codegen('run_quad_hover','-config',cfg,'-args',{D.C(2).kind,D.C(2).G,2,pv,pv,60,S.Ts,E,S.sig,nd},'-o','run_quad_hover_mex');
a=run_quad_hover(D.C(2).kind,D.C(2).G,2,pv,pv,60,S.Ts,E,S.sig,nd);b=run_quad_hover_mex(D.C(2).kind,D.C(2).G,2,pv,pv,60,S.Ts,E,S.sig,nd);assert(isequaln(a.x,b.x)&&isequaln(a.u,b.u));Q=[Q;table("Hover",max(abs(a.x-b.x),[],'all'),'VariableNames',{'Case','MaxStateDifference'})];writetable(Q,fullfile(P.results,'extra_mex_equivalence.csv'));fprintf('EXTRA_MEX_EQUIVALENCE_PASS\n');
end
