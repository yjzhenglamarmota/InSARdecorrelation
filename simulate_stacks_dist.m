% This program aims to simulate redundant interferograms stacking vs independent interferogram
% stacking. Distributed scatterers are assumed.
clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.5 0.5])


rhos=0.8;% correlation between nearest SLCs
snr=10000;numlook=200;numSLC=4;
numsimu=length(rhos);
anglestdind=zeros(numsimu,1);anglestdindnoatm=zeros(numsimu,1);
anglestdred=zeros(numsimu,1);anglestdrednoatm=zeros(numsimu,1);
cmplxstdind=zeros(numsimu,1);cmplxstdindnoatm=zeros(numsimu,1);
cmplxstdred=zeros(numsimu,1);cmplxstdrednoatm=zeros(numsimu,1);

atmphase1=4/3*pi*rand(numSLC,1);    % with atmospheric noise
atmphase=zeros(numSLC,1);           % without atmospheric noise
for k=1:numsimu
%     if mod(k,10)==1
%         disp(k)
%     end
%% SLC simulation
% generate numSLC SLCs, numSLC/2 before and numSLC/2 after. The after scenes will have a
% sysmatic phase shift of 0 rad from the first 10 SLCs (to reduce unwrapping problem)
% note that the phase unwrapping problem cannot be eliminated
rho=rhos(k);
numpixel=2500*numlook;

[SLCs,a]=simulateSLC(numSLC,numpixel,rho,snr);
SLCs1=SLCs.*exp(1j*atmphase1); % SLCs with atm noise
SLCs2=SLCs.*exp(1j*atmphase);  % SLCs without atm noise
% closurephase_comp(numlook,SLCs,1,2,3);

%% create interferograms
for i=1:numSLC/2
    for j=numSLC/2+1:numSLC
        igramsl1=SLCs1(i,:).*conj(SLCs1(j,:));
        igramsl2=SLCs2(i,:).*conj(SLCs2(j,:));
        signalsl=abs(a(i,:)).^2.*rho^(j-i).*exp(1j*(atmphase1(i)-atmphase1(j)));
        noisesl=igramsl1-signalsl;
        signalml=cpxlooks(signalsl,numlook);
        noiseml=cpxlooks(noisesl,numlook);
        igram1{i,j-numSLC/2}=cpxlooks(igramsl1,numlook);
        igram2{i,j-numSLC/2}=cpxlooks(igramsl2,numlook);
        signal{i,j-numSLC/2}=signalml;
        noise{i,j-numSLC/2}=noiseml;
    end
end


%% independent stack
phasestackatm=zeros(1,numpixel/numlook);
phasestacknoatm=zeros(1,numpixel/numlook);
signalstack=zeros(1,numpixel/numlook);noisestack=zeros(1,numpixel/numlook);
anglephasestackatm=zeros(1,numpixel/numlook);
anglephasestacknoatm=zeros(1,numpixel/numlook);
for i=1:numSLC/2
    phasestackatm=phasestackatm+(igram1{i,i});
    phasestacknoatm=phasestacknoatm+(igram2{i,i});
    signalstack=signalstack+signal{i,i};noisestack=noisestack+noise{i,i};
    anglephasestackatm=anglephasestackatm+angle(igram1{i,i});
    anglephasestacknoatm=anglephasestacknoatm+angle(igram2{i,i});
end
anglephasestackatm=anglephasestackatm/numSLC*2;
phasestackatm=angle(phasestackatm);
anglephasestacknoatm=anglephasestacknoatm/numSLC*2;
phasestacknoatm=angle(phasestacknoatm);
%% redundant stack
phasestack2atm=zeros(1,numpixel/numlook);
phasestack2noatm=zeros(1,numpixel/numlook);
signalstack2=zeros(1,numpixel/numlook);noisestack2=zeros(1,numpixel/numlook);
anglephasestack2atm=zeros(1,numpixel/numlook);
anglephasestack2noatm=zeros(1,numpixel/numlook);
for i=1:numSLC/2
    for j=1:numSLC/2
        phasestack2atm=phasestack2atm+(igram1{i,j});
        phasestack2noatm=phasestack2noatm+(igram2{i,j});
        signalstack2=signalstack2+signal{i,j};noisestack2=noisestack2+noise{i,j};
        anglephasestack2atm=anglephasestack2atm+angle(igram1{i,j});
        anglephasestack2noatm=anglephasestack2noatm+angle(igram2{i,j});
    end
end
anglephasestack2atm=anglephasestack2atm/numSLC/numSLC*4;
phasestack2atm=angle(phasestack2atm);
anglephasestack2noatm=anglephasestack2noatm/numSLC/numSLC*4;
phasestack2noatm=angle(phasestack2noatm);

%% save the results
anglestdind(k)=std(anglephasestackatm(:));anglestdindnoatm(k)=std(anglephasestacknoatm(:));
anglestdred(k)=std(anglephasestack2atm(:));anglestdrednoatm(k)=std(anglephasestack2noatm(:));
cmplxstdind(k)=std(phasestackatm(:));cmplxstdindnoatm(k)=std(phasestacknoatm(:));
cmplxstdred(k)=std(phasestack2atm(:));cmplxstdrednoatm(k)=std(phasestack2noatm(:));

%% compute covariance matrice between noise/signal in different interferograms
% independent
% for i=1:numSLC/2
%     noiseind{i}=noise{i,i};
%     signalind{i}=signal{i,i};
% end
% noisecovind=zeros(numSLC/2,numSLC/2);signalcovind=zeros(numSLC/2,numSLC/2);
% noiseind=reshape(noiseind,numSLC/2,1);signalind=reshape(signalind,numSLC/2,1);
% for i=1:numSLC/2
%     for j=1:numSLC/2
%         noisecovind(i,j)=mycoh(noiseind{i},noiseind{j},2500);
%         signalcovind(i,j)=mycoh(signalind{i},signalind{j},2500);
%     end
% end
% figure;imagesc(noisecovind);colorbar;axis image;title('noise covariance, ind')
% figure;imagesc(signalcovind);colorbar;axis image;title('atm covariance, ind')

% % redundant
% noisecov=zeros(numSLC^2/4,numSLC^2/4);signalcov=zeros(numSLC^2/4,numSLC^2/4);
% noise=reshape(noise,numSLC^2/4,1);signal=reshape(signal,numSLC^2/4,1);
% for i=1:numSLC^2/4
%     for j=1:numSLC^2/4
%         noisecov(i,j)=mycoh(noise{i},noise{j},2500);
%         signalcov(i,j)=mycoh(signal{i},signal{j},2500);
%     end
% end
% figure;imagesc(noisecov);colorbar;axis image;title('noise covariance, red')
% figure;imagesc(signalcov);colorbar;axis image;title('signal covariance, red')


end


%% plots
figure;hold on;grid on;
plot(rhos.^(numSLC/2),anglestdind,'k*-','linewidth',2);
plot(rhos.^(numSLC/2),anglestdred,'r*-','linewidth',2);
plot(rhos.^(numSLC/2),anglestdindnoatm,'k*--','linewidth',2);
plot(rhos.^(numSLC/2),anglestdrednoatm,'r*--','linewidth',2);
hold off;
legend('independent (with atm noise)','redundant (with atm noise)','independent (no atm noise)','redundant (no atm noise)')
title('averaging angles')
xlabel('correlation of median temporal-span interferogram');ylabel('std phase')
% saveas(gcf,'avgangle','png')
figure;hold on;grid on;
plot(rhos.^(numSLC/2),cmplxstdind,'k*-','linewidth',2);
plot(rhos.^(numSLC/2),cmplxstdred,'r*-','linewidth',2);
plot(rhos.^(numSLC/2),cmplxstdindnoatm,'k*--','linewidth',2);
plot(rhos.^(numSLC/2),cmplxstdrednoatm,'r*--','linewidth',2);
hold off;
legend('independent (with atm noise)','redundant (with atm noise)','independent (no atm noise)','redundant (no atm noise)')
title('averaging complex numbers')
xlabel('correlation of median temporal-span interferogram');ylabel('std phase')
% saveas(gcf,'avgcomplx','png')