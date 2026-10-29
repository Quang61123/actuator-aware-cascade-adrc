function [Ms,rho,wgrid,Si,wpeak] = chan_ms(kind,gj,Ts,gin,act,bs,oact,nfreq)
% Refined numerical input-sensitivity peak, not a certified continuous bound.
if nargin<6,bs=[];end
if nargin<7,oact=[];end
if nargin<8,nfreq=600;end
M=chan_model(kind,gj,Ts,gin,act,bs,oact);k=M.keep;
A=M.A(k,k);B=M.B(k,7);C=M.Ct(k);D=M.Dt(7);
rho=max(abs(eig(A)));wpeak=NaN;
if rho>=1,Ms=Inf;wgrid=[];Si=[];return,end
sys=ss(A,B,C,D,Ts);wgrid=[0 logspace(-2,log10(pi/Ts),nfreq)];
Si=reshape(freqresp(sys,wgrid),1,[]);mag=abs(Si);[Ms,j]=max(mag);wpeak=wgrid(j);
peaks=find(mag(2:end-1)>=mag(1:end-2)&mag(2:end-1)>=mag(3:end))+1;
peaks=peaks(mag(peaks)>1+1e-8); % Ignore round-off peaks on the unity high-frequency plateau.
opts=optimset('TolX',1e-7,'Display','off');
for j=peaks
 lo=max(wgrid(j-1),1e-9);hi=wgrid(j+1);
 [lw,neg]=fminbnd(@(lw)-abs(evalfr(sys,exp(1i*exp(lw)*Ts))),log(lo),log(hi),opts);
 if -neg>Ms,Ms=-neg;wpeak=exp(lw);end
end
end
