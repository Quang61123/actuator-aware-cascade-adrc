function [co,ledger]=cascade_outer_design(co,S,dub,cap)
original=co;ratios=logspace(log10(1.5),log10(60),20);if co.kind==2,ratios=1;end
candidates=cell(numel(ratios)+1,1);cost=inf(numel(ratios)+1,1);ledger=table();candidates{end}=original;z=cascade_stats(original,S);if valid(original,true,S,dub,cap),cost(end)=z(1);end
parfor j=1:numel(ratios)
 cc=original;if cc.kind==2,shape=original.G(1:5)';upper=4;else,shape=ratios(j);upper=S.wo_max/shape;end
 grid=logspace(log10(.05),log10(upper),36);best=NaN;idx=0;for i=numel(grid):-1:1,cc.G(1:6)=gains_of(cc.kind,shape,grid(i))';cc.G(7:12)=cc.G(1:6);if valid(cc,false,S,dub,cap),best=grid(i);idx=i;break,end,end
 if isfinite(best)
  if idx<numel(grid),lo=log(best);hi=log(grid(idx+1));for b=1:14,w=exp((lo+hi)/2);cc.G(1:6)=gains_of(cc.kind,shape,w)';cc.G(7:12)=cc.G(1:6);if valid(cc,false,S,dub,cap),lo=log(w);else,hi=log(w);end,end;best=exp(lo);end
  cc.G(1:6)=gains_of(cc.kind,shape,best)';cc.G(7:12)=cc.G(1:6);if valid(cc,true,S,dub,cap),zz=cascade_stats(cc,S);cost(j)=zz(1);candidates{j}=cc;end
 end
end
[m,j]=min(cost);assert(isfinite(m),'No feasible cascade design for %s at cap %.8g',co.name,cap);co=candidates{j};
for j=1:numel(candidates),if ~isempty(candidates{j}),zz=cascade_stats(candidates{j},S);ledger=[ledger;table(string(co.name),cap,j,zz(1),zz(4),isfinite(cost(j)),'VariableNames',{'Ctrl','CascadeCap','Candidate','SigmaE','CascadeSigmaDU','Feasible'})];end,end
fprintf('CASCADE_DESIGN %s cap=%.8g sigmae=%.8g sigmaDU=%.8g\n',co.name,cap,m,cascade_stats(co,S)*[0;0;0;1;0]);
end
function ok=valid(c,full,S,dub,cap)
  zz=cascade_stats(c,S);ok=all(isfinite(zz))&&zz(2)<=S.ub(1)*(1+1e-9)&&zz(3)<=dub(1)*(1+1e-9)&&zz(4)<=cap*(1-0.01);if ~ok,return,end
  ms=chan_ms(c.kind,c.G(1:6)',S.Ts,c.G(19:24)',S.act,[],[c.tauhat c.ndhat],600);ok=ms<=S.Msmax;if ~ok,return,end
  if full,os=chan_step(c.kind,c.G(1:6)',S.Ts,c.G(19:24)',S.act,4*S.Tstep(1),[c.tauhat c.ndhat]);else,os=step_screen(c.kind,c.G(1:6)',S.Ts,c.G(19:24)',S.act,4*S.Tstep(1),[c.tauhat c.ndhat]);end;ok=os<=S.osmax;
 end
