% This program plots coh-T
clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
addpath realdata_stat
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','inches')
set(groot, 'defaultFigurePosition',[0 0 10 5.45])

% read igram lists
igramnames=importdata('Cascadia_ulist');
file1= 'redallstack2';file2= 'indallstack2';
% read stacks file of Cascadia
fid=fopen(file1);nr=1104;dat=fread(fid,[nr*2,inf],'float','ieee-le');
amp=dat(1:nr,:)';phaseed=dat(nr+1:end,:)';naz=size(amp,1);fclose(fid);
fid=fopen(file2);nr=1104;dat=fread(fid,[nr*2,inf],'float','ieee-le');
phaseind=dat(nr+1:end,:)';naz=size(amp,1);fclose(fid);
amp(amp==0)=nan;pic=plotamp(amp,nr,naz,15,85);

cohall=[];vphall=[];count=0;posindex=[];
stat_igramlist=struct;
for i=1:length(igramnames)
    igram=igramnames{i};
    daterg=igram(1:17);date1=daterg(1:8);date2=daterg(10:17);
    timespan=datenum(date2,'yyyymmdd')-datenum(date1,'yyyymmdd');
    statname=['stat_' daterg];
    statdata=importdata(statname);    
    % track the same resolution cell over different igrams
    cohphase=statdata.data(:,3);nresol=length(cohphase);
    varphase=statdata.data(:,5);avgphase=statdata.data(:,4);
    for j=1:nresol
        if statdata.data(j,2)<420
            continue;
        end
        loc=[statdata.data(j,1) statdata.data(j,2)];
        locstr=['loc' num2str(loc(1)) num2str(loc(2))];
        if isfield(stat_igramlist, locstr)
            stat_igramlist.(locstr).coh=[stat_igramlist.(locstr).coh; cohphase(j)];
            stat_igramlist.(locstr).timespan=[stat_igramlist.(locstr).timespan;timespan];
            stat_igramlist.(locstr).varphase=[stat_igramlist.(locstr).varphase;varphase(j)];
            stat_igramlist.(locstr).avgphase=[stat_igramlist.(locstr).avgphase;avgphase(j)];
            stat_igramlist.(locstr).igramname=[stat_igramlist.(locstr).igramname;igram];
        else
            stat_igramlist.(locstr)=struct;
            stat_igramlist.(locstr).coh=[];
            stat_igramlist.(locstr).timespan=[];
            stat_igramlist.(locstr).varphase=[];
            stat_igramlist.(locstr).avgphase=[];
            stat_igramlist.(locstr).igramname=[];
            count=count+1;
            posindex=[posindex;count loc(1) loc(2)];
%             figure(1);scatter(loc(1),loc(2),100,'r','filled')
            % store all coherence and phase variance data for later analysis
             cohall=[cohall;cohphase];vphall=[vphall;varphase];
        end
    end         
end
fields=fieldnames(stat_igramlist);
%%%% Yey plots and analysis
%% First, var(\phi) and coherence
numlook=150;
% cohtheo=linspace(0.01,max(cohall)+0.1,1000);
% prevarigram=(1-cohtheo.^2)./2./cohtheo.^2/numlook; % theory
% vphallhi=prctile(vphall,95);vphall(vphall>vphallhi)=nan; % trim down the real data a bit
% figure(2);hold on;scatter(cohall,vphall,'filled');
% % scatter(cohtheo,prevarigram);
% % legend('Real data','theory')
% xlabel('Coherence');ylabel('var(\phi_{decor})');grid on;
% xlim([0.02,1])

%% Next coherence in different igrams
% for i=1:1:length(fieldnames(stat_igramlist))
%     gama=stat_igramlist.(fields{i}).coh;
%     if mean(gama)>0.1
%     T=stat_igramlist.(fields{i}).timespan;
%     [Tsort,index]=sort(T);gama=gama(index);
%     figure(3);subplot(1,2,1);scatter(Tsort,gama,'filled');
%     ylim([0.05,0.6]);xlabel('Time span, days'); ylabel('coherence');grid on;
%     legend(['nr=' num2str(posindex(i,2)) ' naz=' num2str(posindex(i,3))]);
%     subplot(1,2,2);
%     imagesc(pic);axis image;axis off;hold on;scatter(posindex(i,2),posindex(i,3),100,'r','filled')
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     saveas(gcf,['nr=' num2str(posindex(i,2)) '_naz=' num2str(posindex(i,3)) '_coh'],'epsc')
%     end
% end

% for i=1:5:length(fieldnames(stat_igramlist))
%     gama=stat_igramlist.(fields{i}).coh;
%     if mean(gama)<0.1
%     T=stat_igramlist.(fields{i}).timespan;
%     figure(4);subplot(1,2,1);scatter(T,gama,'filled');
%     ylim([0,0.2]);xlabel('Time span, days'); ylabel('coherence');grid on;
%     legend(['nr=' num2str(posindex(i,2)) ' naz=' num2str(posindex(i,3))]);
%     subplot(1,2,2);
%     imagesc(pic);axis image;hold on;axis off;scatter(posindex(i,2),posindex(i,3),100,'r','filled')
%     saveas(gcf,['nr=' num2str(posindex(i,2)) '_naz=' num2str(posindex(i,3)) '_coh'],'png')
%     end
% end

% plot a few slected points for the paper
for i=1:1:length(fieldnames(stat_igramlist))
    gama=stat_igramlist.(fields{i}).coh;
    if posindex(i,2)==506 && posindex(i,3)==1526
        T=stat_igramlist.(fields{i}).timespan;
        [Tsort,index]=sort(T);gama=gama(index);
        figure;scatter(Tsort,gama,'filled');
        ylim([0.05,1]);grid on;
        fig=gcf;
        fig.InvertHardcopy = 'off';
        set(gcf,'color','white')
        saveas(gcf,['CAS_nr=' num2str(posindex(i,2)) '_naz=' num2str(posindex(i,3)) '_coh'],'epsc')
    end
    if posindex(i,2)==146 && posindex(i,3)==2126
        T=stat_igramlist.(fields{i}).timespan;
        [Tsort,index]=sort(T);gama=gama(index);
        figure;scatter(Tsort,gama,'filled');
        ylim([0.05,1]);grid on;
        fig=gcf;
        fig.InvertHardcopy = 'off';
        set(gcf,'color','white')
        saveas(gcf,['CAS_nr=' num2str(posindex(i,2)) '_naz=' num2str(posindex(i,3)) '_coh'],'epsc')
    end
end

%% brief investigation into periodic coherence
% testloc=stat_igramlist.loc6262306;
% goodcohigramlist=testloc.igramname(testloc.coh>0.22 & testloc.timespan>200,:);
% 
lon0=-124.55555725;lat0=48.98333359;
dlon=0.34722223e-4*80; % take looks into account
dlat=-0.13888889e-03*20;

rloc=506;azloc=1526;
lon=lon0+dlon*(rloc-1)
lat=lat0+dlat*(azloc-1)

rloc=146;azloc=2126;
lon=lon0+dlon*(rloc-1)
lat=lat0+dlat*(azloc-1)

