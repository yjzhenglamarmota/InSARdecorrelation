function [COV_DECOR,flag]=mydecorcov(coh_sar,A,varphase,gammainf)
% model based on Piyush and Simons paper
% coh_sar : SAR decorrelation phase corerlation matrix
% A: incidence matrix
% varphase: phase variance of the interforograms as listed in A_ind


%% first, inspect coh_sar and eliminates rows/columns that are zeros
nslc=size(coh_sar,1);
nanloc=isnan(coh_sar);nanlines=(sum(nanloc,1)==nslc-1); % SLCs that simply does not cover this area
coh_sar(nanlines,:)=[];coh_sar(:,nanlines)=[];
A(:,nanlines)=[]; A=A(any(A,2),:);
effigram=(sum(A,2)==0);A=A(effigram,:);
numigram=size(A,1);
numrealigram=length(varphase);
if numigram~=numrealigram
    flag=1;COV_DECOR=[];
    return
end
%% Now we can compute
coh_ifg_piyush=0.5*A*coh_sar*A';
diagcoh=diag(coh_ifg_piyush);
coh_ifg_d=diag(diagcoh);
coh_ifg_d_me=1-(-coh_ifg_d+1).^2;
coh_ifg=2*gammainf*(coh_ifg_piyush-coh_ifg_d)+coh_ifg_d_me;
diagcohifg=diag(coh_ifg);
D_diag=sqrt(varphase)./sqrt(diagcohifg);


D=diag(D_diag);
COV_DECOR=D*coh_ifg*D;

flag=0;
return