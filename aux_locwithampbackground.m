% make amplitude plots with locations
clear; clc; close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
addpath realdata_stat
%% amplitude file
file= 'redallstack2';
% file='stackred_20170330';
fid=fopen(file);nr=1104;dat=fread(fid,[nr*2,inf],'float','ieee-le');
amp=dat(1:nr,:)';naz=size(amp,1);fclose(fid);
% amp=amp(:,1:535);nr=size(amp,2);naz=size(amp,1);
amp(amp==0)=nan;pic=plotamp(amp,nr,naz,15,85);

% point location
locnr=[506;146];locnaz=[1526;2126];
% locnr=[506;326];locnaz=[506;566];
% plot
imagesc(pic);axis image;hold on;axis off;
scatter(locnr,locnaz,200,'rh','filled')

fig=gcf;
fig.InvertHardcopy = 'off';
% saveas(gcf,'CAS_pointloc','epsc')
saveas(gcf,'DV_pointloc','epsc')
% xticks([1 nr])
% xticklabels({num2str(lon0),num2str((nr-1)*dlon+lon0)})
% yticks([-1 -0.8 -0.2 0 0.2 0.8 1])