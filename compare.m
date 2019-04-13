clear;clc;close all;
% compare mean and variance between two stacks
addpath /Users/yjzheng/Documents/MATLAB/mytools/
addpath /Users/yjzheng/Documents/MATLAB/bigsentinelscenetest/
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.3 0.3])
file1= 'ssestacked_20160120';
file2= 'ssestackind_20160120';
nr=1104;
%% read data
[amp1,phase1]=readunwfile(file1,nr);
[amp2,phase2]=readunwfile(file2,nr);
naz=size(phase1,1);
addpath ccfiles
[amp3,cc]=readunwfile('ccstack_sse1_ed',nr);
mask=zeros(naz,nr);
amp1(amp1==0)=nan;
% amplow=prctile(amp1(:),20);amphi=prctile(amp1(:),85);amp1(isnan(amp1))=0;
% amp1(amp1<amplow)=amplow;amp1(amp1>amphi)=amphi;
% ampscale=(amp1-amplow)./(amphi-amplow);
mask(isnan(amp1))=1;
phase1(mask==1)=nan;
phase2(mask==1)=nan;
% figure;imagesc(phase1);
%% choose data points 
xdata=50:50:nr-50;
ydata=50:50:naz-50;
[x,y]=meshgrid(xdata,ydata);
xsize=size(x,1);
ysize=size(x,2);

k=0;
expphase1 =zeros(xsize*ysize,1);
expphase2 =zeros(xsize*ysize,1);
varphase1 =zeros(xsize*ysize,1);
varphase2 =zeros(xsize*ysize,1);
ccout     =zeros(xsize*ysize,1);
for i=1:xsize
    for j=1:ysize
        k=k+1;
    	phaseout1=phase1(y(i,j)-20:y(i,j)+20,x(i,j)-20:x(i,j)+20);
        phaseout2=phase2(y(i,j)-20:y(i,j)+20,x(i,j)-20:x(i,j)+20);
        phaseout1=phaseout1(:);phaseout2=phaseout2(:);
        ccoutall=cc(y(i,j)-10:y(i,j)+10,x(i,j)-10:x(i,j)+10);
        ccout(k)=nanmean(ccoutall(:));
        expphase1(k)=nanmean(phaseout1);expphase2(k)=nanmean(phaseout2);
        varphase1(k)=nanvar(phaseout1);varphase2(k)=nanvar(phaseout2);
    end
end
[ccout,I]=sort(ccout);
expphase1=expphase1(I);expphase2=expphase2(I);
varphase1=varphase1(I);varphase2=varphase2(I);
%% make plots
figure;hold on;scatter(ccout,expphase2,'filled');scatter(ccout,expphase1,'filled');
hold off;legend('independent','redundant');xlim([0.02,0.4]);ylim([-5,5]);grid on;
xlabel('Coherence');ylabel('Mean phase');saveas(gcf,'compphasemean','png')

figure;hold on;scatter(ccout,varphase2,'filled');scatter(ccout,varphase1,'filled');
hold off;legend('independent','redundant');xlim([0.02,0.4]);ylim([0,2]);grid on;
xlabel('Coherence');ylabel('Phase variance');saveas(gcf,'compphasevar','png')