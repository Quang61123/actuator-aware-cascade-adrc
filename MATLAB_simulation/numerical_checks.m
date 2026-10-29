function numerical_checks()
P=sim_config();D=load(fullfile(P.results,'design_gains.mat'));C=D.C;S=D.S;pv=quad_params();H=table();
for c=1:4
 E=make_exo(6,91001,40,S.Ts,S);pc=quad_params([],C(c).tauhat,C(c).ndhat);a=run_quad_sim(C(c).kind,C(c).G,6,pv,pc,40,S.Ts,E,S.sig,1);b=run_quad_sim(C(c).kind,C(c).G,6,pv,pc,40,S.Ts,E,S.sig,2);assert(all(isfinite(b.x),'all'));aa=run_metrics(a,[10 40],pv);bb=run_metrics(b,[10 40],pv);delta=max(abs(a.x-b.x),[],'all');H=[H;table(string(C(c).name),delta,(bb.rmse3-aa.rmse3)*1000,'VariableNames',{'Ctrl','MaxStateDifference','RMSEChange_mm'})];
 if c==2,q=run_quad_stress(C(c).kind,C(c).G,6,pv,pc,40,S.Ts,E,S.sig,5*ones(40000,1));assert(isequal(a.x,q.x)&&isequal(a.wcmd,q.wcmd),'History-buffer reduction failed');save(fullfile(P.results,'nominal_waveform.mat'),'a','C');end
end
writetable(H,fullfile(P.results,'rk4_check.csv'));save(fullfile(P.results,'rk4_check.mat'),'H');
% Independent augmented-state expression for increment covariance in all nominal classes.
A=table();ch=[1 3 4 6];for c=1:4,for h=1:4
 g=C(c).G(6*(ch(h)-1)+(1:6))';gi=[];syi=0;svi=0;if h==1,gi=C(c).G(19:24)';syi=S.sy(3);svi=S.sv(3);end;oa=[C(c).tauhat C(c).ndhat];M=chan_model(C(c).kind,g,S.Ts,gi,S.act,[],oa);k=M.keep;F=M.A(k,k);Bn=M.B(k,2:5);cu=M.Cu(k);du=M.Du(2:5);W=diag([S.sy(h)^2 S.sv(h)^2 syi^2 svi^2]);Pn=dlyap(F,Bn*W*Bn');Fa=[F zeros(numel(k),1);cu 0];Ba=[Bn;du];Ca=[cu -1];Pa=dlyap(Fa,Ba*W*Ba');v=Ca*Pa*Ca'+du*W*du';[~,~,~,ref]=chan_var(C(c).kind,g,S.Ts,S.sy(h),S.sv(h),0,gi,syi,svi,S.act,oa);err=abs(v-ref)/max(ref,eps);assert(err<1e-8);res=norm(Pn-F*Pn*F'-Bn*W*Bn','fro')/max(norm(Pn,'fro'),eps);A=[A;table(string(C(c).name),h,err,res,'VariableNames',{'Ctrl','Class','IncrementRelativeError','LyapunovRelativeResidual'})];
end,end
writetable(A,fullfile(P.results,'covariance_check.csv'));
% Fresh nonlinear noise-only hover validation, paired seeds and no bias/wind.
created=false;if isempty(gcp('nocreate')),parpool('Processes',P.workers);created=true;end;guard=onCleanup(@()close_owned(created));C=[C D.Strict];indices=[1 2 4 8 9 10];R=zeros(60,10);
parfor job=1:60
 c=indices(mod(job-1,6)+1);q=floor((job-1)/6)+1;E=make_exo(2,96000+q,60,S.Ts,S);E.bias(:)=0;E.wind(:)=0;if exist('run_quad_hover_mex','file')==3,hover=@run_quad_hover_mex;else,hover=@run_quad_hover;end;L=hover(C(c).kind,C(c).G,2,pv,quad_params([],C(c).tauhat,C(c).ndhat),60,S.Ts,E,S.sig,5*ones(60000,1));assert(all(isfinite(L.x),'all')&&~any(L.wcmd>=pv(11)-1e-6|L.wcmd<=1e-6,'all'));v=L.u(L.t>=20,:);ii=[1 3 4 6];R(job,:)=[var(v(:,ii),1,1),var(diff(v(:,ii)),1,1),var(v(:,5),1),var(diff(v(:,5)),1)];
end
R=reshape(R,[6 10 10]);N=table();K=table();
for j=1:6,c=indices(j);for h=1:4
 g=C(c).G(6*(ch(h)-1)+(1:6))';gi=[];syi=0;svi=0;if h==1,gi=C(c).G(19:24)';syi=S.sy(3);svi=S.sv(3);end;oa=[C(c).tauhat C(c).ndhat];[~,~,vu,vdu]=chan_var(C(c).kind,g,S.Ts,S.sy(h),S.sv(h),0,gi,syi,svi,S.act,oa);a=reshape(R(j,:,h),10,1);d=reshape(R(j,:,h+4),10,1);N=[N;table(string(C(c).name),h,mean(a)/vu,mean(d)/vdu,sqrt(mean(a)),sqrt(mean(d)),'VariableNames',{'Ctrl','Class','VarianceRatioU','VarianceRatioDU','MeasuredSigmaU','MeasuredSigmaDU'})];
end
 M=chan_model_inner(C(c).kind,C(c).G(1:6)',S.Ts,C(c).G(19:24)',S.act,[],oa);Q=chan_model(C(c).kind,C(c).G(1:6)',S.Ts,C(c).G(19:24)',S.act,[],oa);assert(isequal(M.A,Q.A)&&isequal(M.B,Q.B));k=M.keep;F=M.A(k,k);Bn=M.B(k,2:5);cu=M.Ci(k);du=M.Di(2:5);W=diag([S.sy(1)^2 S.sv(1)^2 S.sy(3)^2 S.sv(3)^2]);Pn=dlyap(F,Bn*W*Bn');vu=cu*Pn*cu'+du*W*du';h=cu*(F-eye(numel(k)));qq=cu*Bn-du;vdu=h*Pn*h'+qq*W*qq'+du*W*du';a=reshape(R(j,:,3),10,1);d=reshape(R(j,:,7),10,1);K=[K;table(string(C(c).name),sqrt(vdu),sqrt(mean(d)),mean(a)/vu,mean(d)/vdu,D.dub(3),'VariableNames',{'Ctrl','PredSigmaDU','MeasuredSigmaDU','CascadeRatioU','CascadeRatioDU','LocalDULimit'})];
end
K.CascadeCap=[repmat(D.cap,3,1);repmat(D.strict,3,1)];K.PitchMeasuredSigmaDU=zeros(6,1);for j=1:6,K.PitchMeasuredSigmaDU(j)=sqrt(mean(R(j,:,10)));end;K.MeasuredPass=max(K.MeasuredSigmaDU,K.PitchMeasuredSigmaDU)<=K.CascadeCap*(1+1e-9);assert(all(K.MeasuredPass(2:end)),'Measured hover cascade cap exceeded');
HoverRaw=table();for j=1:6,for seed=1:10,HoverRaw=[HoverRaw;table(string(C(indices(j)).name),96000+seed,sqrt(R(j,seed,7)),sqrt(R(j,seed,10)),K.CascadeCap(j),'VariableNames',{'Ctrl','Seed','MeasuredSigmaDU','PitchMeasuredSigmaDU','CascadeCap'})];end,end;writetable(HoverRaw,fullfile(P.results,'hover_per_seed.csv'));
save(fullfile(P.results,'noise_only.mat'),'R','N','K','C','S');writetable(N,fullfile(P.results,'noise_only.csv'));writetable(K,fullfile(P.results,'cascade_noise.csv'));fprintf('NUMERICAL_CHECKS_PASS\n');
end
function close_owned(created),if created,p=gcp('nocreate');if ~isempty(p),delete(p);end,end,end
