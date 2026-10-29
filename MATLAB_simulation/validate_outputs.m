function validate_outputs()
P=sim_config();D=load(fullfile(P.results,'trials.mat'));U=readtable(fullfile(P.results,'summary.csv'),'TextType','string');assert(height(D.T)==1830&&sum(U.Runs)==1830);assert(sum(U.Failures)==sum(D.T.Failed));assert(sum(U.ThresholdExceedances)==sum(D.T.FluctuationThreshold));assert(all(D.T.Delay<=20)&&all(D.T.Tau>=.015&D.T.Tau<=.12));assert(all(U.RMSE_mean_mm(U.Failures<U.Runs)>=0));
G=load(fullfile(P.results,'design_gains.mat'));assert(numel(G.C)==7&&numel(G.Strict)==3);assert(~any(contains(string({G.C.name}),'-R')));for j=1:numel(G.Curve),assert(isfinite(G.Curve(j).Gamma));end
if exist('run_quad_stress_mex','file')==3
X=table();cases=[5 .03 5;2 .12 10;4 .12 10];for j=1:size(cases,1)
 c=cases(j,1);pv=quad_params([],cases(j,2),cases(j,3));pc=quad_params([],G.C(c).tauhat,G.C(c).ndhat);E=make_exo(6,92001,40,G.S.Ts,G.S);nd=cases(j,3)*ones(40000,1);a=run_quad_stress(G.C(c).kind,G.C(c).G,6,pv,pc,40,G.S.Ts,E,G.S.sig,nd);b=run_quad_stress_mex(G.C(c).kind,G.C(c).G,6,pv,pc,40,G.S.Ts,E,G.S.sig,nd);assert(isequal(isfinite(a.x),isfinite(b.x)));ok=isfinite(a.x)&isfinite(b.x);err=max(abs(a.x(ok)-b.x(ok)));assert(err<1e-6);X=[X;table(string(G.C(c).name),cases(j,2)*1000,cases(j,3),err,'VariableNames',{'Ctrl','Tau_ms','Delay_ms','MaxFiniteStateDifference'})];
end
writetable(X,fullfile(P.results,'mex_off_nominal_equivalence.csv'));
else
 warning('sim:NoMEX','MEX unavailable: using MATLAB simulations; live MEX equivalence check skipped.');
end
CA=readtable(fullfile(P.results,'cascade_design_audit.csv'));assert(height(CA)==10&&all(CA.PredSigmaDU<=CA.CascadeCap));HV=readtable(fullfile(P.results,'hover_per_seed.csv'),'TextType','string');assert(height(HV)==60);US=U(U.Experiment=="S5",:);assert(height(US)==4&&sum(US.Runs)==120);assert(sum(U.Runs(U.Experiment=="E7"))==720);
CV=readtable(fullfile(P.results,'covariance_check.csv'));assert(height(CV)==16&&all(CV.IncrementRelativeError<1e-8));SL=readtable(fullfile(P.results,'simulink.csv'));assert(height(SL)==4&&all(SL.MaxPositionDiff<1e-8));M=readtable(fullfile(P.results,'mex_equivalence.csv'));assert(height(M)==4&&all(M.MaxStateDifference<1e-6));RK=readtable(fullfile(P.results,'rk4_check.csv'));assert(height(RK)==4);
files={'run_pipeline.m','sim_config.m','primary_design.m','mc_trials.m','mc_summarize.m','mc_paired.m','linear_checks.m','numerical_checks.m','validate_outputs.m','make_figures.m','build_mex.m','ctrl_law_logged.m','run_quad_hover.m','chan_model_inner.m'};for k=1:numel(files),q=checkcode(files{k},'-id');for j=1:numel(q),assert(~contains(q(j).id,'PARSE'));end,end
fprintf('VALIDATION_PASS: 1830 raw/summary records; counts; 16 independent covariance cases; 4 MEX/RK4/Simulink cases; source syntax\n');
end
