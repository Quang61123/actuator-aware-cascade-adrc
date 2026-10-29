function q=channel_stats(kind,gj,c,S,gin,oact,plants,fullStep,nfreq)
if nargin<8,fullStep=false;end
if nargin<9,nfreq=600;end
if isempty(gin),[vd,vn,vu,du]=chan_var(kind,gj,S.Ts,S.sy(c),S.sv(c),S.qd(c),[],0,0,S.act,oact);
else,[vd,vn,vu,du]=chan_var(kind,gj,S.Ts,S.sy(c),S.sv(c),S.qd(c),gin,S.sy(3),S.sv(3),S.act,oact);end
q=[sqrt(vd+vn),sqrt(vu),sqrt(du),0,0,0];
if ~all(isfinite(q(1:3))),q(4:6)=Inf;return,end
for p=1:size(plants,1)
 [ms,rho]=chan_ms(kind,gj,S.Ts,gin,plants(p,:),[],oact,nfreq);
 q(4)=max(q(4),ms);q(6)=max(q(6),rho);
end
if fullStep,q(5)=chan_step(kind,gj,S.Ts,gin,S.act,4*S.Tstep(c),oact);
else,q(5)=step_screen(kind,gj,S.Ts,gin,S.act,4*S.Tstep(c),oact);end
end
