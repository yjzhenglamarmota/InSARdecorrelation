% This program aims to simulate redundant interferograms stacking vs independent interferogram
% stacking. Distributed scatterers are assumed.
clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.5 0.5])


rhos=0.9:-0.05:0.5;% correlation between nearest SLCs
snr=1000000;numlook=200;numSLC=10;
numsimu=length(rhos);
atmphase1=4/3*pi*rand(numSLC,1);    % with atmospheric noise
atmphase=zeros(numSLC,1);           % without atmospheric noise
count=0;
for k =1:numsimu
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

% closurephase_comp(numlook,SLCs,1,2,3);

%% create interferograms
for i=1:numSLC/2
    for j=numSLC/2+1:numSLC        
        igramsl=SLCs(i,:).*conj(SLCs(j,:));     
        igram{i,j-numSLC/2}=cpxlooks(igramsl,numlook);        
    end
end

%% Let's find out the relation beween phi_decor vs coherence
% figure;
for i=1:numSLC/2
    for j=numSLC/2+1:numSLC
        count=count+1;
        gama(count)=mycoh(SLCs(i,:),SLCs(j,:),2500);
        igramnow=angle(igram{i,j-numSLC/2});
        stdigram(count)=std(igramnow(:));
%         subplot(2,2,count);histogram(igramnow(:));title(['std= ' num2str(stdigram(count))])
    end
end
if k==1
% figure;igramnow=angle(igram{5,1});histogram(igramnow(:));xlim([-pi,pi])
% title(['mean=' num2str(mean(igramnow(:))) 'std=' num2str(std(igramnow(:)))]);
% gama1=mycoh(SLCs(5,:),SLCs(6,:),numlook);
% saveas(gcf,['igramcoh= ' num2str(gama1)],'png');
% figure;igramnow=angle(igram{3,1});histogram(igramnow(:));xlim([-pi,pi])
% title(['mean=' num2str(mean(igramnow(:))) 'std=' num2str(std(igramnow(:)))]);
% gama2=mycoh(SLCs(3,:),SLCs(6,:),numlook);
% saveas(gcf,['igramcoh= ' num2str(gama2)],'png');

% figure;igramnow=angle(igram{1,1});histogram(igramnow(:));xlim([-pi,pi])
% title(['mean=' num2str(mean(igramnow(:))) 'std=' num2str(std(igramnow(:)))]);
% gama3=mycoh(SLCs(1,:),SLCs(6,:),numlook);
% saveas(gcf,['igramcoh= ' num2str(gama3)],'png');

% noisecoh=zeros(numSLC^2/4,numSLC^2/4);
% noise={igram{1,1},igram{2,2},igram{2,1},igram{1,2}};
% for i=1:numSLC^2/4
%     for j=1:numSLC^2/4
%         noisecoh(i,j)=mycoh(angle(noise{i}),angle(noise{j}),2500);
%     end
% end
% figure;imagesc(noisecoh);colorbar;axis image;title('noise coherence matrix')
% cohstackind=[0.5 0.5]*noisecoh(1:2,1:2)*[0.5 0.5]';
% cohstackred=[0.25 0.25 0.25 0.25]*noisecoh*[0.25 0.25 0.25 0.25]';
% D=[stdigram(1), stdigram(4), ...
%     stdigram(3), stdigram(2)];
% D=diag(D);
% noisecov=D*noisecoh*D;
% figure;imagesc(noisecov);colorbar;axis image;title('noise covariance matrix')
% covstackind=[0.5 0.5]*noisecov(1:2,1:2)*[0.5 0.5]';
% covstackred=[0.25 0.25 0.25 0.25]*noisecov*[0.25 0.25 0.25 0.25]';

end
end

figure;scatter(gama,stdigram.^2,'filled');grid on;
xlabel('Coherence');ylabel('var(\phi_{decor})')
xlim([0.02,0.4])
hold on;
prevarigram=(1-gama.^2)./2./gama.^2/numlook;
scatter(gama,prevarigram);
legend('simulation','theory')
% save decorphasedata gama varigram