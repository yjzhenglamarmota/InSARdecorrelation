% analyze two phasetriplets
clear; close all;clc;
addpath /Users/yjzheng/Documents/MATLAB/mytools/

[amp,phase_sl]=readflatfile('igram_134_sl_looks',1200);
[~,phase_ml]=readflatfile('igram_134_ml',1200);

% figure
% subplot(1,3,1)
% imagesc(amp);axis image;axis off;colorbar;
% amplo=prctile(amp(:),10);amphi=prctile(amp(:),90);
% caxis([amplo,amphi])
% subplot(1,3,2)
% imagesc(phase_sl);axis image;colorbar;axis off;
% subplot(1,3,3)
% imagesc(phase_ml);axis image;colorbar;axis off;

% figure
% [pic,phlow,phhi]=plotphase(phase_ml,amp,1200,size(amp,2),5,95,1,99);
% imagesc(pic);axis image;axis off;colorbar;caxis([phlow,phhi]);colormap jet;

%% mountain,farm and city
mountain=phase_ml(71:120,140:189);
% figure
% imagesc((mountain));colorbar;axis off;axis image;
% title('mountain')

farm=phase_ml(585:585+49,809:809+49);
% figure
% imagesc((farm));colorbar;axis off;axis image;
% title('farm')

city=phase_ml(987:987+49,1024:1024+49);
% figure
% imagesc((city));colorbar;axis off;axis image;
% title('city')

%% histograms
figure;histogram(mountain(:));
title(['mountain mean=' num2str(mean(mountain(:))) ' std=' num2str(std(mountain(:)))]);
saveas(gcf,'mountain_hist','png')
figure;histogram(farm(:));
title(['farm mean=' num2str(mean(farm(:))) ' std=' num2str(std(farm(:)))]);
saveas(gcf,'farm_hist','png')
figure;histogram(city(:));
title(['city mean=' num2str(mean(city(:))) ' std=' num2str(std(city(:)))]);
saveas(gcf,'city_hist','png')