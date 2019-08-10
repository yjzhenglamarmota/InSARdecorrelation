% This program simulates a stack of coherence and varphase
% and compares predictions from three decorrelation covariance model
clear;clc;close all;
set(0,'defaultAxesFontSize', 25);
set(groot, 'defaultFigureUnits','inches')
set(groot, 'defaultFigurePosition',[0 0 10 10])
timespan=3*365; % observation period, days
rvtime=12;      % satellite revisit time
numslc=floor(timespan/rvtime); % number of acquisitions
numindstack=floor(numslc/2);   % number of igrams in an independent stack
numredstack=numindstack.^2;    % number of igrams in a redundant stack
numslc=numindstack*2;

%% construct A matrices
I1=eye(numindstack);I2=-eye(numindstack);
numperm=0;
for i=1:numindstack
    if numperm<1
        A=[I1 I2];
    else
        I3=[I2(numperm+1:end,:);I2(1:numperm,:)];
        A=[A;I1 I3];
    end
    numperm=numperm+1;
end
A_ind=A(1:numindstack,:);
A_red=A;

deltats=zeros(numredstack,1);
for i=1:numredstack
    Aline=A(i,:);
    indexloc=find(Aline);
    deltats(i)=(indexloc(2)-indexloc(1))*rvtime;
end
%%
gamainf=0.1; 
taulist=1:1:100;
% taulist=1:1000:10001;
ind_hanssen=zeros(length(taulist),1);red_hanssen=zeros(length(taulist),1);
ind_agram=zeros(length(taulist),1); red_agram=zeros(length(taulist),1);
ind_yujie=zeros(length(taulist),1); red_yujie=zeros(length(taulist),1); redindratioyujie=zeros(length(taulist),1);
for iter=1:length(taulist)
%     gamainf=gamainflist(iter);tau=taulong;
    tau=taulist(iter);
%% construct coh_sar matrix
coh_sar=zeros(numslc,numslc);
for i=1:numslc
    for j=i:numslc
        dt=(j-i)*rvtime;
        [coh_sar(i,j)]=simulatecoh(gamainf,tau,dt);
        coh_sar(j,i)=coh_sar(i,j);
    end
end

%% construct varphase
[cohall,varphaseall]=simulatecoh(gamainf,tau,deltats);
% varphaseall=ones(size(varphaseall));
varphaseind=varphaseall(1:numindstack);
varphasered=varphaseall;

%% visualize the coherence-t plot
% coh_plot=simulatecoh(gamainf,tau,0:10:1200);
% varphase_plot=(1-coh_plot.^2)./2./coh_plot.^2;
% 
% figure(2);hold on; plot(0:10:1200,varphase_plot,'linewidth',3);
% xlabel('Time span, days');ylabel('Phase variance \sigma_{\phi}^2, [rad]^2');grid on;
% 
% figure(3);hold on; plot(0:10:1200,coh_plot,'linewidth',3)
% hold off;xlabel('Time span, days'); ylabel('Coherence \rho');grid on;
% 
%% Predictions of three models
% Hanssen (2001)
covdecorind_hanssen=diag(varphaseind); varpreind_hanssen=mean(covdecorind_hanssen(:));
covdecorred_hanssen=diag(varphasered); varprered_hanssen=mean(covdecorred_hanssen(:));
% Agram and Simons (2015)
covdecorind_agram=piyushdecorcov(coh_sar,A_ind,varphaseind);varpreind_agram=mean(covdecorind_agram(:));
covdecorred_agram=piyushdecorcov(coh_sar,A_red,varphasered);varprered_agram=mean(covdecorred_agram(:));
% Yujie (2019)
covdecorind_yujie=mydecorcov(coh_sar,A_ind,varphaseind,gamainf);varpreind_yujie=mean(covdecorind_yujie(:));
covdecorred_yujie=mydecorcov(coh_sar,A_red,varphasered,gamainf);varprered_yujie=mean(covdecorred_yujie(:));

ind_hanssen(iter)=varpreind_hanssen;
red_hanssen(iter)=varprered_hanssen;
ind_agram(iter)=varpreind_agram;
red_agram(iter)=varprered_agram;
ind_yujie(iter)=varpreind_yujie;
red_yujie(iter)=varprered_yujie;
redindratioyujie(iter)=varprered_yujie/varpreind_yujie;
end

% figure(1); hold on;
% dhanssen=plot(taulist,ind_hanssen,'g','linewidth',3);
% dagram=plot(taulist,ind_agram,'b','linewidth',3);
% dyujie=plot(taulist,ind_yujie,'r','linewidth',3);
% plot(taulist,red_hanssen,'g--','linewidth',2);
% plot(taulist,red_agram,'b--','linewidth',2);
% plot(taulist,red_yujie,'r--','linewidth',2)
% indeg=plot(taulist,taulist*nan,'k','linewidth',3);
% redeg=plot(taulist,taulist*nan,'k--','linewidth',2);
% % xlabel('\rho_\infty')
% xlabel('\tau, days')
% ylabel('Predicted phase noise \sigma_{\phi}^2, [rad]^2')
% grid on
% a=axes('position',get(gca,'position'),'visible','off');
% legend([dhanssen,dagram,dyujie],'Hanssen', 'Agram-Simons','Proposed')
% legend(a,[indeg,redeg],'Independent','Redundant')
% fig=gcf;fig.InvertHardcopy = 'off';
% set(gcf,'color','white')
% saveas(gcf,'varydecor_compare3models_test2','epsc')

% figure(2)
% 
% legend('\tau=1 day','\tau=100 days','\tau=1000 days')
% fig=gcf;fig.InvertHardcopy = 'off';
% set(gcf,'color','white')
% saveas(gcf,'varphi_t_wtvarytau','epsc')
% % 
% figure(3)
% legend('\tau=1 day','\tau=100 days','\tau=1000 days')
% fig=gcf;fig.InvertHardcopy = 'off';
% set(gcf,'color','white')
% saveas(gcf,'rho_t_wtvarytau','epsc')

figure(4)
scatter(taulist,redindratioyujie)

function [gama,varphase]=simulatecoh(gamainf,tau, deltats)
% generate a list of coherence and phase variance
% gama=gamainf+(1-gamainf)*exp(-t/tau)
% input: gamainf,tau: parameters of the coherecne-t function
% input: deltats: list of time spans
% output: gama: list of coherence corresponding to the list of time spans
% output: varphase: list of phase variance w.r.t. the list of time spans

gama=gamainf+(1-gamainf).*exp(-deltats./tau);
varphase=(1-gama.^2)./2./gama.^2;

return
end