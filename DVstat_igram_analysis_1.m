% Derived from stat_igram_analysis_1.m
% For death valley region
clear;clc;close all;
addpath /Users/yjzheng/iCloud/Documents/MATLAB/mytools/
addpath realdata_stat2
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','inches')
set(groot, 'defaultFigurePosition',[0 0 10 10])

%% read igram lists
igramnames=importdata('DV_ulist');
redigramnames=importdata('sselisted_20170330');
indigramnames=importdata('sselistind_20170330');

%% make slc list, initialize coh_sar matrix
nslc=2*length(indigramnames);
coh_sar=nan(nslc,nslc);
for i=1:nslc
    coh_sar(i,i)=1;
end
slcnames=cell(nslc,1);
for i=1:length(indigramnames)
    igram=indigramnames{i};
    daterg=igram(1:17);date1=daterg(1:8);date2=daterg(10:17);
    slcnames{i}=date1; slcnames{i+nslc/2}=date2;
end
%% read stacks file of death valley
nr=720;
file1= 'stackred_20170330';file2= 'stackind_20170330';
fid=fopen(file1);dat=fread(fid,[nr*2,inf],'float','ieee-le');
amp=dat(1:nr,:)';phaseed=dat(nr+1:end,:)';naz=size(amp,1);fclose(fid);
fid=fopen(file2);dat=fread(fid,[nr*2,inf],'float','ieee-le');
phaseind=dat(nr+1:end,:)';naz=size(amp,1);fclose(fid);
amp(amp==0)=nan;pic=plotamp(amp,nr,naz,15,85);
% figure('Name','1','position',[0,0,0.5,0.5]);hold on; 
% imagesc(pic);set(gca,'Ydir','reverse');axis image;
cohall=[];vphall=[];count=0;posindex=[];
stat_igramlist=struct;
%% read stat_igram files and write structure stat_igramlist
for i=1:length(igramnames)
    igram=igramnames{i};
    redflag=sum(strcmp(igram,redigramnames));
    indflag=sum(strcmp(igram,indigramnames));
    daterg=igram(1:17);date1=daterg(1:8);date2=daterg(10:17);
    rows=strcmp(date1,slcnames); rowindx=find(rows==1);
    cols=strcmp(date2,slcnames); colindx=find(cols==1);
    if sum(rows)<1 || sum(cols)<1
        continue;
    end
    timespan=datenum(date2,'yyyymmdd')-datenum(date1,'yyyymmdd');
    statname=['stat_' daterg];
    statdata=importdata(statname);    
    % track the same resolution cell over different igrams
    cohphase=statdata.data(:,3);nresol=length(cohphase);
    varphase=statdata.data(:,5);avgphase=statdata.data(:,4);
    cohall=[cohall;cohphase];vphall=[vphall;varphase];
    for j=1:nresol
        if (statdata.data(j,2)<360 && statdata.data(j,1)>360) % ignore the upper right part of the igram (wrong DEM)
            continue;
        end
        loc=[statdata.data(j,1) statdata.data(j,2)];
        locstr=['loc' num2str(loc(1)) '_' num2str(loc(2))];
        if isfield(stat_igramlist, locstr)
            stat_igramlist.(locstr).coh=[stat_igramlist.(locstr).coh; cohphase(j)];
            stat_igramlist.(locstr).timespan=[stat_igramlist.(locstr).timespan;timespan];
            stat_igramlist.(locstr).varphase=[stat_igramlist.(locstr).varphase;varphase(j)];
            stat_igramlist.(locstr).avgphase=[stat_igramlist.(locstr).avgphase;avgphase(j)];
            stat_igramlist.(locstr).neffflag=[stat_igramlist.(locstr).neffflag;indflag];
            stat_igramlist.(locstr).redflag=[stat_igramlist.(locstr).redflag;redflag];
            stat_igramlist.(locstr).cohsar(rowindx,colindx)=cohphase(j);
            stat_igramlist.(locstr).cohsar(colindx,rowindx)=cohphase(j);
            stat_igramlist.(locstr).igram{end+1}=igram;
        else
            stat_igramlist.(locstr)=struct;
            stat_igramlist.(locstr).coh=cohphase(j);
            stat_igramlist.(locstr).timespan=timespan;
            stat_igramlist.(locstr).varphase=varphase(j);
            stat_igramlist.(locstr).avgphase=avgphase(j);
            stat_igramlist.(locstr).neffflag=indflag;
            stat_igramlist.(locstr).redflag=redflag;
            stat_igramlist.(locstr).cohsar=coh_sar;
            stat_igramlist.(locstr).cohsar(rowindx,colindx)=cohphase(j);
            stat_igramlist.(locstr).cohsar(colindx,rowindx)=cohphase(j);
            stat_igramlist.(locstr).igram{1}=igram;
            count=count+1;
            posindex=[posindex;count loc(1) loc(2)];
%             figure(1);scatter(loc(1),loc(2),100,'r','filled')          
        end
    end         
end
fields=fieldnames(stat_igramlist);

vphallhi=prctile(vphall,90);
disp('Files read and loaded.')
%% construct A matrix
A_ind=zeros(length(indigramnames),nslc);A_red=zeros(length(redigramnames),nslc);
for i=1:length(indigramnames)
    igram=indigramnames{i};
    daterg=igram(1:17);date1=daterg(1:8);date2=daterg(10:17);
    rows=strcmp(date1,slcnames); rowindx=find(rows==1);
    cols=strcmp(date2,slcnames); colindx=find(cols==1);
    A_ind(i,rowindx)=1;A_ind(i,colindx)=-1;
end
for i=1:length(redigramnames)
    igram=redigramnames{i};
    daterg=igram(1:17);date1=daterg(1:8);date2=daterg(10:17);
    rows=strcmp(date1,slcnames); rowindx=find(rows==1);
    cols=strcmp(date2,slcnames); colindx=find(cols==1);
    A_red(i,rowindx)=1;A_red(i,colindx)=-1;
end
%% compare phase variance between stacks and individual igrams
winsize=25;neffind=zeros(count,1);neffred=zeros(count,1);
gamainf=zeros(count,1);
prered_piyush=nan(count,1); preind_piyush=nan(count,1);
prered=nan(count,1); preind=nan(count,1);
prered_hanssen=nan(count,1); preind_hanssen=nan(count,1);
varphasered=nan(count,1); varphaseind=nan(count,1);
for i=1:count
    if mod(i,10)==0, disp(i); end
    %% stacks variance from data
        x=posindex(i,2);y=posindex(i,3);
        phaseout1= phaseed(y-winsize+1:y+winsize,x-winsize+1:x+winsize);
        phaseout2= phaseind(y-winsize+1:y+winsize,x-winsize+1:x+winsize);
        phaseout1=phaseout1(:);phaseout2=phaseout2(:);
        expphase1=nanmean(phaseout1);expphase2=nanmean(phaseout2);
        varphase1=nanvar(phaseout1);varphase2=nanvar(phaseout2);
    %% phase variance of individual igrams
        varphaseigram=stat_igramlist.(fields{i}).varphase;
        varphaseigram(varphaseigram>vphallhi)=nan;
        locigramnames=stat_igramlist.(fields{i}).igram;
        gama=stat_igramlist.(fields{i}).coh;
        tspn=stat_igramlist.(fields{i}).timespan;
        % rough estimate first
        param=pinv([ones(length(gama),1),tspn])*log(gama);
        lnA=param(1);invtao=param(2);
        f= @(b,x) b(1)+(1-b(1)).*exp(b(2).*x);
        B=fminsearchbnd(@(b) norm(gama-f(b,tspn)),[1-exp(lnA),-1/invtao],[0,-1],[1,0]);
%          figure(2);scatter(tspn,gama,'pg');hold on;scatter(tspn,f(B,tspn),'filled');ylim([0,1]);
%          close(figure(2));
        gamainf(i)=B(1);
        indflag=logical(stat_igramlist.(fields{i}).neffflag);
        neffind(i)=sum(indflag);
        redflag=logical(stat_igramlist.(fields{i}).redflag);
        neffred(i)=sum(redflag);
        varphaseigramind=varphaseigram(indflag);
        neffind(i)=nnz(isfinite(varphaseigramind));
        varphasemeanind=nansum(varphaseigramind)/neffind(i);
        varphaseigramred=varphaseigram(redflag);
        neffred(i)=nnz(isfinite(varphaseigramred));
        varphasemeanred=nansum(varphaseigramred)/neffred(i);       
        
     %% construct decorrelation covariance matrix
        coh_sar=stat_igramlist.(fields{i}).cohsar;
        % piyush model for independent stack
        COV_DECOR_piyush_ind=piyushdecorcov(coh_sar,A_ind,varphaseigramind);
        COV_DECOR_piyush_red=piyushdecorcov(coh_sar,A_red,varphaseigramred);
        COV_DECOR_ind=mydecorcov(coh_sar,A_ind,varphaseigramind,gamainf(i));
        COV_DECOR_red=mydecorcov(coh_sar,A_red,varphaseigramred,gamainf(i));
   
        M=neffind(i);
        if M>=20
            
            % redundant stack
            prered_piyush(i)=nanmean(COV_DECOR_piyush_red(:));
            prered(i)=nanmean(COV_DECOR_red(:));
            prered_hanssen(i)=varphasemeanred/neffred(i);
            % independent           
            preind_piyush(i)=nanmean(COV_DECOR_piyush_ind(:));
            preind(i)=nanmean(COV_DECOR_ind(:));
            preind_hanssen(i)=varphasemeanind/M;
          
            varphaseind(i)=varphase2;
            varphasered(i)=varphase1;

        end
end
%% plot(s)
figure(8); hold on;
sz=140;
scatter(varphasered,varphaseind,sz,'k','linewidth',2,'MarkerFaceColor',rgb('beige')); % data
scatter(prered,preind,sz+20,'p','filled',... % yujie model
    'MarkerEdgeColor',[0.5, 0, 0.5],'MarkerFaceColor',[0.7,0,0.7],'linewidth',2);
scatter(prered_piyush, preind_piyush, sz, 'd',... % piyush model
    'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0,0.7,0.7],'linewidth',2);
scatter(prered_hanssen,preind_hanssen,sz,'s','linewidth',2,... % hanssen model
   'MarkerEdgeColor',rgb('orange'),'MarkerFaceColor',rgb('orange'));
% add y=x line to plots
yy=0:0.1:0.6;xx=yy;
plot(xx,yy,'r--','linewidth',2);
% xlabel('\sigma^2(\phi_{red}^{decor}),[rad]^2');ylabel('\sigma^2(\phi_{ind}^{decor}),[rad]^2');
grid on; axis image;
xlim([0,0.6]);ylim([0,0.6]);
[~, objh] = legend( {'Data','Proposed','Agram-Simons','Hanssen'},'location','Southeast','Fontsize',25);
objhl= findobj(objh, 'type', 'patch'); %// children of legend of type line
set(objhl, 'Markersize', 25); %// set value as desired
fig=gcf;fig.InvertHardcopy = 'off';set(gcf,'color','white')
saveas(gcf,'DV_validate','epsc')
