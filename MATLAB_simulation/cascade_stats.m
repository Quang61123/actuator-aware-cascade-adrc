function z=cascade_stats(co,S)
M=chan_model_inner(co.kind,co.G(1:6)',S.Ts,co.G(19:24)',S.act,[],[co.tauhat co.ndhat]);k=M.keep;F=M.A(k,k);B=M.B(k,2:5);W=diag([S.sy(1)^2 S.sv(1)^2 S.sy(3)^2 S.sv(3)^2]);rho=max(abs(eig(F)));z=inf(1,5);z(5)=rho;if rho>=1,return,end
P=dlyap(F,B*W*B');c=M.Cu(k);d=M.Du(2:5);h=c*(F-eye(numel(k)));q=c*B-d;vu=c*P*c'+d*W*d';vdu=h*P*h'+q*W*q'+d*W*d';ce=M.Cy(k);bd=M.A(k,M.id);xss=(eye(numel(k))-F)\bd;if abs(ce*xss)>1e-7*max(1,norm(ce)*norm(xss)),return,end
vd=S.qd(1)*S.Ts*(xss'*dlyap(F',ce'*ce)*xss);ve=ce*P*ce'+vd;ci=M.Ci(k);di=M.Di(2:5);hi=ci*(F-eye(numel(k)));qi=ci*B-di;vi=hi*P*hi'+qi*W*qi'+di*W*di';z=[sqrt(max(0,ve)),sqrt(max(0,vu)),sqrt(max(0,vdu)),sqrt(max(0,vi)),rho];
end
