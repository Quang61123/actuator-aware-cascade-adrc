function primary_design()
P=sim_config();D=load(fullfile(P.results,'initial_designs.mat'));S=D.S;C=D.C;dub=D.dub;Curve=D.Curve;K=readtable(fullfile(P.results,'initial_cascade_noise.csv'),'TextType','string');cap=K.MeasuredSigmaDU(K.Ctrl=="PID-Budget");strict=dub(3);A=table();
for c=[2 4 5 6 7]
 dl=dub;if c==7,dl=inf(1,4);end;[C(c),L]=cascade_outer_design(C(c),S,dl,cap);A=[A;L];
end
save(fullfile(P.results,'primary_design.mat'),'C','S','dub','Curve','cap','strict','A');
Strict=C([1 2 4]);for j=1:3
 Strict(j).name=[Strict(j).name '-Strict'];if j==1,sc=.9;for ch=[3 4 5 6],ii=6*(ch-1)+(1:6);g=Strict(j).G(ii);g(1:4)=g(1:4).*[sc^2;sc^3;sc;sc];Strict(j).G(ii)=g;end;else,sc=.8;for ch=[3 4 5 6],ii=6*(ch-1)+(1:6);Strict(j).G(ii(1:2))=sc*Strict(j).G(ii(1:2));end;end
 [Strict(j),L]=cascade_outer_design(Strict(j),S,dub,strict);A=[A;L];
end
for j=1:numel(Curve),[Curve(j).Controller,L]=cascade_outer_design(Curve(j).Controller,S,Curve(j).Gamma*dub,cap);A=[A;L];end
save(fullfile(P.results,'design_gains.mat'),'C','S','dub','Curve','Strict','cap','strict','P');save(fullfile(P.results,'design_search.mat'),'A');writetable(A,fullfile(P.results,'cascade_search.csv'));finalize_design();fprintf('CASCADE_DESIGN_PASS\n');
end
