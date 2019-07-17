function [COV_DECOR,flag]=mydecorcov_test(coh_sar,A,varphase,gammainf)
% model described in Yujie's Redundant interferogram stacking paper
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
coh_ifg=zeros(size(A,1));
for i=1:size(A,1)
    Aline=A(i,:);
    ind1=find(Aline);
    for j=i:size(A,1)
        Aline=A(j,:);
        ind2=find(Aline);
        rho_13=coh_sar(ind1(1),ind2(1));
        rho_24=coh_sar(ind1(2),ind2(2));
        coh_ifgtemp=rho_13*rho_24-gammainf^2;
        coh_ifgtemp=coh_ifgtemp/(1-gammainf^2);
        coh_ifg(i,j)=1-(1-coh_ifgtemp).^0.45;
        coh_ifg(j,i)=coh_ifg(i,j);
        
    end
end


D_diag=sqrt(varphase);
D=diag(D_diag);
COV_DECOR=D*coh_ifg*D;

flag=0;
return