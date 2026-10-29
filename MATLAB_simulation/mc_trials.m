function mc_trials()
P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));C=D.C;S=D.S;B=D.Curve;Jobs=struct('Experiment',{},'Scenario',{},'Ctrl',{},'Controller',{},'Seed',{},'Tau',{},'Delay',{},'Noise',{},'Pattern',{},'Gamma',{});
for c=1:7,ex='E1';if c>4,ex='E3';end;for q=1:30,add(ex,C(c),P.nominalSeed+q,.03,5,1,'fixed',NaN);end,end
for b=1:numel(B),for q=1:30,add('E4',B(b).Controller,P.nominalSeed+q,.03,5,1,'fixed',B(b).Gamma);end,end
for c=1:4,for z=1:6,for q=1:10
 scale=[1 .5 2 4 1 1];pattern={'fixed','fixed','fixed','fixed','sine','burst'};add(['E5-' num2str(z)],C(c),P.stressSeed+q,.03,5,scale(z),pattern{z},NaN);
end,end,end
for c=[1 2 4],for tau=P.taus,for nd=P.delays,for q=1:10,add('E7',C(c),P.gridSeed+q,tau,nd,1,'fixed',NaN);end,end,end,end
for c=1:numel(D.Strict),for q=1:30,add('E0',D.Strict(c),P.nominalSeed+q,.03,5,1,'fixed',NaN);end,end
for c=1:4,for q=1:30,add('S5',C(c),P.nominalSeed+q,.03,5,1,'fixed',NaN);end,end
save(fullfile(P.results,'job_protocol.mat'),'Jobs','P');useMex=exist('run_quad_stress_mex','file')==3;fields={'rmse3','chat_cmd','rotor_std','sat_rotor','sat_virt','max3'};batchSize=30;nr=numel(Jobs);chunks=ceil(nr/batchSize);R=nan(nr,7);
for chunk=1:chunks
 fn=fullfile(P.results,sprintf('batch_%03d.mat',chunk));ii=(chunk-1)*batchSize+1:min(chunk*batchSize,nr);
 if isfile(fn),A=load(fn);assert(isequaln(A.J,Jobs(ii)),'Checkpoint protocol mismatch');R(ii,:)=A.Z;continue,end
 J=Jobs(ii);Z=nan(numel(ii),7);
 parfor k=1:numel(ii)
  j=J(k);pv=quad_params([],j.Tau,j.Delay);pc=quad_params([],j.Controller.tauhat,j.Controller.ndhat);duration=40;win=[10 40];if j.Scenario==5,duration=45;win=[5 45];end;E=make_exo(j.Scenario,j.Seed,duration,S.Ts,S);n=round(duration/S.Ts);t=(0:n-1)'*S.Ts;nd=j.Delay*ones(n,1);if strcmp(j.Pattern,'sine'),nd=round(5+5*sin(2*pi*t));elseif strcmp(j.Pattern,'burst'),nd(t>=10&mod(t-10,2)<.2)=20;end
  if useMex&&j.Scenario~=5,L=run_quad_stress_mex(j.Controller.kind,j.Controller.G,j.Scenario,pv,pc,duration,S.Ts,E,S.sig*j.Noise,nd);elseif j.Scenario==5&&exist('run_quad_s5_mex','file')==3,L=run_quad_s5_mex(j.Controller.kind,j.Controller.G,j.Scenario,pv,pc,duration,S.Ts,E,S.sig*j.Noise,nd);else,L=run_quad_stress(j.Controller.kind,j.Controller.G,j.Scenario,pv,pc,duration,S.Ts,E,S.sig*j.Noise,nd);end;v=run_metrics(L,win,pv);Z(k,:)=[v.rmse3*1000,v.chat_cmd,v.rotor_std,v.sat_rotor*100,v.sat_virt*100,v.diverged,v.max3*1000];
 end
 tmp=[fn '.partial'];save(tmp,'J','Z','fields','-v7');movefile(tmp,fn,'f');R(ii,:)=Z;fprintf('MC %d/%d completed\n',max(ii),nr);
end
T=struct2table(rmfield(Jobs,'Controller'));T.RMSE_mm=R(:,1);T.Increment=R(:,2);T.RotorStd=R(:,3);T.RotorClip_pct=R(:,4);T.VirtualClip_pct=R(:,5);T.Failed=R(:,6)>0;T.MaxError_mm=R(:,7);T.FluctuationThreshold=R(:,3)>20&~T.Failed;
save(fullfile(P.results,'trials.mat'),'Jobs','R','T','P','C','S','fields');writetable(T,fullfile(P.results,'trials.csv'));mc_summarize();mc_paired();fprintf('TRIALS_PASS: %d saved realizations\n',nr);
 function add(ex,co,seed,tau,nd,noise,pat,gamma)
 Jobs(end+1)=struct('Experiment',ex,'Scenario',6-strcmp(ex,'S5'),'Ctrl',co.name,'Controller',co,'Seed',seed,'Tau',tau,'Delay',nd,'Noise',noise,'Pattern',pat,'Gamma',gamma);
 end
end
