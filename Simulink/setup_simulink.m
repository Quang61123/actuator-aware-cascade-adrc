function setup_simulink(scen,controller,duration)
if nargin<1,scen=6;end;if nargin<2,controller='AA-Budget';end;if nargin<3,duration=40;end
here=fileparts(mfilename('fullpath'));root=fullfile(here,'..','MATLAB_simulation');addpath(root,here);P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));C=D.C;S=D.S;i=find(strcmp({C.name},controller),1);assert(~isempty(i));Ts=S.Ts;N=round(duration/Ts);pvt=quad_params();if scen==4,pvt=quad_params([1.3 1.3 2],1.5*pvt(16),pvt(17));end;pv=quad_params([],C(i).tauhat,C(i).ndhat);E=make_exo(scen,P.nominalSeed+1,duration,Ts,S);
values=struct('Ts',Ts,'Tend',duration-Ts,'scen',scen,'kind',C(i).kind,'G',C(i).G,'pv',pv,'pvt',pvt,'sig',S.sig,'exo_ts',[(0:N-1)'*Ts,[E.noise;E.bias;E.wind]']);names=fieldnames(values);for k=1:numel(names),assignin('base',names{k},values.(names{k}));end
end
