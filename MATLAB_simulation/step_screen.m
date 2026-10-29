function os=step_screen(kind,gj,Ts,gin,act,Tend,oact)
% Exact lifting for a fast 10-ms preliminary screen; final gains use chan_step.
M=chan_model(kind,gj,Ts,gin,act,[],oact);k=M.keep;A=M.A(k,k);B=M.B(k,6);C=M.Cy(k);
if max(abs(eig(A)))>=1,os=Inf;return,end
stride=10;As=A^stride;Bs=zeros(size(B));term=B;
for j=1:stride,Bs=Bs+term;term=A*term;end
N=ceil(Tend/(Ts*stride));y=ltitr(As,Bs,ones(N,1))*C';
os=max(0,100*(max(y)-1));
end
