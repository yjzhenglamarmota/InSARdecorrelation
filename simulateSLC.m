function SLCs=simulateSLC(numslc,numpixel,cc)
% model for distributed scatterers (with )
% a_(i+1)=cc*a_i+(1-cc^2)^0.5*d_(i+1),a_1=c
% d is a circular Gaussian random variable ~ N(0,0.5)
% c is a constant complex number. For distributed scatterers, can be
% initiated as a circular Gaussian.
% input: 
% numslc: number of SLCs wanted
% numpixel: number of pixels in each SLC
% cc: a list of correlation with respect to time. Should be numslc in
% length (e.g. [1, e^{t1/tau}, e^{t2/tau},...])

% output:
% SLCs is a complex numslc by numpixel matrix, each row represents one SLC.
% Alternative interpretation, each column represents a different realizatio
% and each row represents a different snapshot in time
mu=0;sigma=sqrt(0.5);
a=zeros(1,numpixel);
SLCs=zeros(numslc,numpixel);

%% create the master SLC
c=normrnd(0,sigma,[1,numpixel*2]);
real=c(:,1:2:end);imag=c(:,2:2:end);c_cp=complex(real,imag);

a(1,:)=c_cp;masterSLC=a(1,:);
SLCs(1,:)=masterSLC;

%% create subsequent SLCs
d        =  normrnd(mu,sigma,[numslc-1,numpixel*2]);
real=d(:,1:2:end);imag=d(:,2:2:end);d_cp=complex(real,imag);

% slave SLCs
cc(cc>1)=1;

for i=2:numslc 
    rho=cc(i)/cc(i-1);
    a(i,:)=rho*a(i-1,:)+sqrt(1-rho^2)*d_cp(i-1,:);
    slaveSLC=a(i,:);
    SLCs(i,:)=slaveSLC;
end
   
return