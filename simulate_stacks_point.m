% This program aims to simulate redundant interferograms stacking vs independent interferogram
% stacking. Point scatterers are assumed.
clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.5 0.5])

%% SLC simulation
% generate 20 SLCs, 10 before and 10 after. The after scenes will have a
% sysmatic phase shift of 3 rad (not sure how much to assign)
numSLC=20;numpixel=250000;
defphase=1.5;
defcmplx=complex(cos(defphase),sin(defphase));
SLCs=simulateSLC(numSLC,numpixel*2,'point',0.02);
SLCs(11:end,:)=SLCs(11:end,:).*defcmplx;

%stat_slcs(SLCs(1:10,:));
%stat_slcs(SLCs(11:end,:));

%% create interferograms
numlook=100;
for i=1:numSLC/2
    for j=11:numSLC
        igramsl=SLCs(i,:).*conj(SLCs(j,:));
        igram{i,j-10}=cpxlooks(igramsl,numlook);
    end
end

%% independent stack
phasestack=zeros(1,numpixel/numlook);
for i=1:numSLC/2
    phasestack=phasestack+angle(igram{i,i});
end
phasestack=phasestack/numSLC*2;
stack=reshape(phasestack,50,50);
figure;imagesc(stack);colorbar;caxis([-pi,pi]);axis image;colormap jet
title(['phase mean=' num2str(mean(phasestack(:))) ' std=' num2str(std(phasestack(:)))]);
saveas(gcf,'ind_stack_point_amp1','png')
%% redundant stack
phasestack2=zeros(1,numpixel/numlook);
for i=1:numSLC/2
    for j=1:numSLC/2
        phasestack2=phasestack2+angle(igram{i,j});
    end
end
phasestack2=phasestack2/numSLC/numSLC*4;
stack=reshape(phasestack2,50,50);
figure;imagesc(stack);colorbar;caxis([-pi,pi]);axis image;colormap jet
title(['phase mean=' num2str(mean(phasestack2(:))) ' std=' num2str(std(phasestack2(:)))]);
saveas(gcf,'rdn_stack_point_amp1','png')

%% stack phase difference
phasediff=phasestack-phasestack2;
stack=reshape(phasestack2,50,50);
figure;imagesc(stack);colorbar;caxis([-pi,pi]);axis image;
title(['phase mean=' num2str(mean(phasediff(:))) ' std=' num2str(std(phasediff(:)))]);
saveas(gcf,'dif_stack_dist','png')
figure;histogram(phasediff(:));
title(['phase mean=' num2str(mean(phasediff(:))) ' std=' num2str(std(phasediff(:)))]);
saveas(gcf,'difhist_stack_dist','png')