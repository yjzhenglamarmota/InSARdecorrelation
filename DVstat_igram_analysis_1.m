% Derived from stat_igram_analysis_1.m
% For death valley region
clear;clc;close all;
addpath /Users/yjzheng/Documents/MATLAB/mytools/
addpath realdata_stat2
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','normalized')
set(groot, 'defaultFigurePosition',[0 0 0.6 0.5])

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
    ohall=[cohall;cohphase];vphall=[vphall;varphase];
    for j=1:nresol
%         if (statdata.data(j,2)<360 && statdata.data(j,1)>360) % ignore the upper right part of the igram (wrong DEM)
%             continue;
%         end
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
for i=1:count
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
%         % rough estimate first
%         param=pinv([ones(length(gama),1),tspn])*log(gama);
%         lnA=param(1);invtao=param(2);
%         f= @(b,x) b(1)+(1-b(1)).*exp(b(2).*x);
%         B=fminsearchbnd(@(b) norm(gama-f(b,tspn)),[0,-1/20],[0,-1/20],[1,0]);
%         figure(2);scatter(tspn,gama,'pg');hold on;scatter(tspn,f(B,tspn),'filled');ylim([0,1]);
%         close(figure(2));
%         gamainf(i)=B(1);
        gamainf(i)=prctile(gama,1);
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
        
%         varphaseigram(varphaseigram>vphallhi)=nan;
%         varphaseigramind=varphaseigram(indflag);
%         varphaseigramred=varphaseigram(redflag);
%         indxigramind=zeros(length(locigramind),1);
%         for j=1:neffind(i)
%             indxigram=strcmp(locigramind(j),locigramred);
%             indxigramind(j)=find(indxigram==1);
%         end
     %% construct decorrelation covariance matrix
        coh_sar=stat_igramlist.(fields{i}).cohsar;
        % piyush model for independent stack
        COV_DECOR_piyush_ind=piyushdecorcov(coh_sar,A_ind,varphaseigramind);
        COV_DECOR_piyush_red=piyushdecorcov(coh_sar,A_red,varphaseigramred);
        COV_DECOR_ind=mydecorcov(coh_sar,A_ind,varphaseigramind,gamainf(i));
        COV_DECOR_red=mydecorcov(coh_sar,A_red,varphaseigramred,gamainf(i));
%         COV_DECOR_ind_test=COV_DECOR_red(indxigramind,:);
%         COV_DECOR_ind_test=COV_DECOR_ind_test(:,indxigramind);
     %% plot
        M=neffind(i);
        if M>=20
            
            % redundant stack
            figure(5);hold on; 
            prered_piyush=nanmean(COV_DECOR_piyush_red(:));
            prered=nanmean(COV_DECOR_red(:));
%             prered=(2*(M-1)*gamainf(i)/(1+gamainf(i))+1)/M/M*varphasemeanred;
            scatter(varphase1,prered_piyush,'b','filled')
            scatter(varphase1,prered,'r','filled')
            scatter(varphase1,varphasemeanred/neffred(i),'pg')
            xlabel('var(\phi_{stack_{red}})');ylabel('predicted var(\phi_{stack_{red}})');
            axis image; xlim([0,1]);ylim([0,1]);grid on;
            
            % independent
            figure(6);hold on;
            preind_piyush=nanmean(COV_DECOR_piyush_ind(:));
             preind=nanmean(COV_DECOR_ind(:));
%               preind=varphasemeanind/M;
            scatter(varphase2,preind_piyush,'b','filled')
            scatter(varphase2,preind,'r','filled')
            scatter(varphase2,varphasemeanind/M,'pg')
            xlabel('var(\phi_{stack_{ind}})');ylabel('predicted var(\phi_{stack_{ind}})');
            axis image;xlim([0,1]);ylim([0,1]);grid on;
       
            figure(7);hold on;
            scatter(varphase1,varphase2,'b','filled');grid on; xlim([0,0.5]);ylim([0,0.5]);
            xlabel('var(\phi_{stack_{red}})');ylabel('var(\phi_{stack_{ind}})');
        
        
            figure(8);hold on;
            prediff_piyush=preind_piyush-prered_piyush;
            prediff=preind-prered;
           
            scatter(varphase2-varphase1,prediff,'r','filled');grid on;
            scatter(varphase2-varphase1,prediff_piyush,'b','filled')
            scatter(varphase2-varphase1,varphasemeanind/neffind(i)-varphasemeanred/neffred(i),'pg')
            xlim([-0.25,0.25]);ylim([-0.25,0.25])
            xlabel('var(\phi_{stack_{diff}})');ylabel('predicted difference');
        
        end
end
% add y=x line to plots
yy=-1:0.1:1;xx=yy;
figure(5);plot(xx,yy,'r','linewidth',2);
figure(6);plot(xx,yy,'r','linewidth',2);
figure(7);plot(xx,yy,'r','linewidth',2);
figure(8);plot(xx,yy,'r','linewidth',2);
