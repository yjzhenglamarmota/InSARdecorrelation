function [SLCs,a]=simulateSLC(numslc,numpixel,cc,snr)
% model:(constant spatial and temporal coherence) for distributed scatterers
% yi=ai+1/sqrt(snr)*nthermal, a_(i+1)=cc*a_i+(1-cc^2)^0.5*d_(i+1),a_1=c
% nthermal,d are circular Gaussian random variables ~ N(0,0.5)
% c is a constant complex number. For distributed scatterers, can be
% initiated as a circular Gaussian.
% input: 
% numslc: number of SLCs wanted
% numpixel: number of pixels in each SLC
% cc: coherence between adjacent a's
% snr: signal to noise ratio between a and nthermal
% The final coherence gama=cc/(1+1/snr)

% output:
% SLCs is a complex numslc by numpixel matrix, each row represents one SLC.

mu=0;sigma=sqrt(0.5);
r=1./sqrt(snr);
a=zeros(1,numpixel);
SLCs=zeros(numslc,numpixel);

%% create the master SLC
c=normrnd(0,sigma,[1,numpixel*2]);
real=c(:,1:2:end);imag=c(:,2:2:end);c_cp=complex(real,imag);
nthermal1=normrnd(mu,sigma,[1,numpixel*2]);
real=nthermal1(:,1:2:end);imag=nthermal1(:,2:2:end);nthermal1_cp=complex(real,imag);
a(1,:)=c_cp;masterSLC=a(1,:)+nthermal1_cp*r;
compphase=-0.01;
compcmplx=complex(cos(compphase),sin(compphase));
% 
% masterSLC(abs(angle(masterSLC)-pi)<0.01)=masterSLC(abs(angle(masterSLC)-pi)<0.01).*compcmplx;
% masterSLC(abs(angle(masterSLC)+pi)<0.01)=masterSLC(abs(angle(masterSLC)+pi)<0.01).*conj(compcmplx);

SLCs(1,:)=masterSLC;

%% create subsequent SLCs
d        =  normrnd(mu,sigma,[numslc-1,numpixel*2]);
nthermal2=  normrnd(mu,sigma,[numslc-1,numpixel*2]);
real=d(:,1:2:end);imag=d(:,2:2:end);d_cp=complex(real,imag);
real=nthermal2(:,1:2:end);imag=nthermal2(:,2:2:end);nthermal2_cp=complex(real,imag);

% slave SLCs
cc(cc>1)=1;
for i=2:numslc 
    a(i,:)=cc*a(i-1,:)+sqrt(1-cc^2)*d_cp(i-1,:);
    slaveSLC=a(i,:)+r.*nthermal2_cp(i-1,:);
%     slaveSLC(abs(angle(slaveSLC)-pi)<0.01)=slaveSLC(abs(angle(slaveSLC)-pi)<0.01).*compcmplx;
%     slaveSLC(abs(angle(slaveSLC)+pi)<0.01)=slaveSLC(abs(angle(slaveSLC)+pi)<0.01).*conj(compcmplx);

    SLCs(i,:)=slaveSLC;
end
   
return