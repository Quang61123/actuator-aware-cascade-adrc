function [g,q,factor]=feasibility_repair(kind,g,c,S,gin,dub,R,nfreq)
% Reduce bandwidth and, for PID only, set-point weight until full verification passes.
if nargin<9,nfreq=4000;end
original=g;factor=1;
for j=1:70
 g=scale_design(kind,original,factor);
 q=channel_stats(kind,g,c,S,gin,R.oact,R.plants,true,nfreq);
 if kind==2 && isfinite(q(5)) && q(5)>S.osmax
  for ib=1:8
   g(5)=g(5)*.8;
   q(5)=chan_step(kind,g,S.Ts,gin,S.act,4*S.Tstep(c),R.oact);
   if q(5)<=S.osmax,break,end
  end
  original(5)=g(5); % Set-point weight changes reference response, not noise/sensitivity.
 end
 if all(isfinite(q))&&q(2)<=S.ub(c)*(1+1e-8)&&q(3)<=dub*(1+1e-8)&&q(4)<=S.Msmax&&q(5)<=S.osmax&&q(6)<1,return,end
 if j<=10,factor=factor*.985;else,factor=factor*.9;end
end
error('No verified feasible gain after bandwidth repair, class %d',c);
end
