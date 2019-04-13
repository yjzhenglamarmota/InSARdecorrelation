clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
addpath /Users/yjzheng/Documents/MATLAB/bigsentinelscenetest/
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.5 0.5])
file1= 'ssestacked_20160120';
file2= 'ssestackind_20160120';
nr=1104;
%% compare plots
% plotSSEstack(file1);saveas(gcf,file1,'png')
% plotSSEstack(file2);saveas(gcf,file2,'png')

%% compare
% [amp1,phase1]=readunwfile(file1,nr);
% [amp2,phase2]=readunwfile(file2,nr);
% naz=size(phase1,1);
% mask=zeros(naz,nr);
% mask(amp1==0)=1;
% phase1(mask==1)=nan;phase1=phase1(420:end,:);
% phase2(mask==1)=nan;phase2=phase2(420:end,:);
% figure;
% subplot(1,2,1);imagesc(phase1);colormap jet; colorbar;caxis([-5,5]);colorbar;axis image;
% subplot(1,2,2);imagesc(phase2);colormap jet; colorbar;caxis([-5,5]);colorbar;axis image;
% saveas(gcf,'stack_compare_phase','png')
% difphase=phase2-phase1;
% amp1=amp1(420:end,:);
% figure;imagesc(difphase);colormap jet; colorbar;caxis([-5,5]);colorbar;axis image;
% [pic,phlow,phhi,mask]=plotphase(difphase,amp1,nr,naz-419,15,85,5,95);
% f1=imagesc(pic);axis image; axis off; colormap jet;colorbar;
% caxis([phlow,phhi]); saveas(gcf,'stackdifphase','png')
% alpha(f1,0.5);[ii,jj]=find(abs(difphase)>4);hold on;scatter(jj,ii,5,'k')
% saveas(gcf,'largestphasemisfitloc1','png')

%% profiles
[amp1,phase1]=readunwfile(file1,nr);
[amp2,phase2]=readunwfile(file2,nr);
naz=size(phase1,1);
mask=zeros(naz,nr);
amp1(amp1==0)=nan;
amplow=prctile(amp1(:),20);amphi=prctile(amp1(:),85);amp1(isnan(amp1))=0;
amp1(amp1<amplow)=amplow;amp1(amp1>amphi)=amphi;
ampscale=(amp1-amplow)./(amphi-amplow);
mask(ampscale==0)=1;
phase1(mask==1)=nan;
phase2(mask==1)=nan;

lines=250:200:850;
linesshift=lines+420;
figure('position',[0,0,0.25,0.75]);hold on;
for i=1:length(lines)
    subplot(length(lines),1,i)
    scatter(1:nr,phase2(linesshift(i),:),20,'filled');hold on;
    scatter(1:nr,phase1(linesshift(i),:),20,'filled');hold off;
    legend('independent','redundant')
    title(['Line' num2str(lines(i))])
end

%% compare with respect to ccfiles
addpath ccfiles
[amp1,cc1]=readunwfile('ccstack_sse1_ed',nr);cc1(cc1==0)=nan;
cc1=cc1(420:end,:);cclo=prctile(cc1(:),5);cchi=prctile(cc1(:),95);
ccgroup=zeros(size(cc1));

ccselect=cc1>0.025;
ccgroup(ccselect)=cc1(ccselect);
ccgroup=sparse(ccgroup);
[ii,jj,value]=find(ccgroup);
difphasegroup1=difphase(ccselect);
figure;
f1=imagesc(cc1);colorbar;axis image;caxis([cclo cchi]);alpha(f1,0.5);
hold on;scatter(jj,ii,5,'r*')
figure
f1=imagesc(cc1);colorbar;axis image;caxis([cclo cchi]);hold on;axis off;
[ii,jj]=find(abs(difphase)>4);f2=scatter(jj,ii,3,'k');alpha(f2,0.3)
saveas(gcf,'largestphasemifitloc2','png')


figure;scatter(value,difphasegroup1);xlabel('cc');ylabel('phase difference, rad');
hold on;
plot(linspace(0,1,100),-4*ones(1,100),'k--','linewidth',3);
plot(linspace(0,1,100),4*ones(1,100),'k--','linewidth',3);

%% 