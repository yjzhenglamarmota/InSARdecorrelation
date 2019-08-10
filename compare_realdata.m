clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
addpath /Users/yjzheng/Documents/MATLAB/bigsentinelscenetest
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.5 0.5])
file1= 'ssestacked_20160120';
file2= 'ssestackind_20160120';
nr=1104;
%% compare plots
file1_masked=[file1 '_masked'];
file2_masked=[file2 '_masked'];

plotSSEstack(file1_masked);saveas(gcf,file1_masked,'epsc')
plotSSEstack(file2_masked);saveas(gcf,file2_masked,'epsc')

%% profiles
% [amp1,phase1]=readunwfile(file1_masked,nr);
% [amp2,phase2]=readunwfile(file2_masked,nr);
% naz=size(phase1,1);
% figure;
% subplot(1,2,1);imagesc(phase1);colormap jet; colorbar;caxis([-5,5]);colorbar;axis image;
% subplot(1,2,2);imagesc(phase2);colormap jet; colorbar;caxis([-5,5]);colorbar;axis image;
% mask=zeros(naz,nr);
% amp1(amp1==0)=nan;
% amplow=prctile(amp1(:),20);amphi=prctile(amp1(:),85);amp1(isnan(amp1))=0;
% amp1(amp1<amplow)=amplow;amp1(amp1>amphi)=amphi;
% ampscale=(amp1-amplow)./(amphi-amplow);
% mask(ampscale==0)=1;
% phase1(mask==1)=nan;
% phase2(mask==1)=nan;

% lines=1150:200:1750;
% figure('position',[0,0,0.4,0.4]);hold on;
% for i=1:length(lines)
% %     subplot(2,2,i)
%     scatter(1:nr,phase2(lines(i),:),20,'filled');hold on;
%     scatter(1:nr,phase1(lines(i),:),20,'filled');hold off;
%     if i==2
%         legend('Independent stacking','Redundant stacking')
%     end
% %     title(['Line' num2str(lines(i))])
%     grid on; ylim([-4,4]);xlim([0,1200])
% %     xlabel('Pixels'); ylabel('Phase, [rad]')
% %     saveas(gcf,['plots/Line' num2str(lines(i))],'epsc')
% end

%% compare phase variance
% lines=1100:100:1800;
% pixels=51:100:nr;
% winsize=50;
% for i=1:length(lines)
%     for j=1:length(pixels)
%         indphase=phase2(lines(i)-winsize/2:lines(i)+winsize/2-1, pixels(j)-winsize/2:pixels(j)+winsize/2-1);
%         redphase=phase1(lines(i)-winsize/2:lines(i)+winsize/2-1, pixels(j)-winsize/2:pixels(j)+winsize/2-1);
%         varphaseind(i,j)=nanvar(indphase(:));
%         varphasered(i,j)=nanvar(redphase(:));
%     end
% end
% varphaseline_ind=mean(nanmean(varphaseind,2));
% varphaseline_red=mean(nanmean(varphasered,2));
% varphaseline_red/varphaseline_ind