function finalize_design()
P=sim_config();fn=fullfile(P.results,'design_gains.mat');D=load(fn);R=struct('plants',D.S.act,'oact',D.S.act);for j=1:numel(D.Strict)
 c=D.Strict(j);for h=2:4,ch=[1 3 4 6];ii=6*(ch(h)-1)+(1:6);g=c.G(ii)';z=channel_stats(c.kind,g,h,D.S,[],R.oact,R.plants,true,P.nfreq);if z(5)>10||z(4)>2||z(2)>D.S.ub(h)||z(3)>D.dub(h),[g,~,~]=feasibility_repair(c.kind,g,h,D.S,[],D.dub(h),R,P.nfreq);c.G(ii)=g';if h==3,c.G(25:30)=g';end;fprintf('STRICT_INNER_REPAIR %s class=%d\n',c.name,h);end;end
 [D.Strict(j),L]=cascade_outer_design(c,D.S,D.dub,D.strict);writetable(L,fullfile(P.results,sprintf('strict_final_search_%d.csv',j)));
end;save(fn,'-struct','D');audit_design();fprintf('CASCADE_DESIGN_PASS\n');
end
