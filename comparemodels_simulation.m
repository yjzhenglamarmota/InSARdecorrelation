% This program simulates a stack of coherence and varphase
% and compares predictions from three decorrelation covariance model
clear;clc;close all;

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
% gamainflst=0.1:0.05:0.95; taushort=10; taulong=600;
gamainf=0.1; taulist=10:50:600;
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
[~,varphaseall]=simulatecoh(gamainf,tau,deltats);
varphaseind=varphaseall(1:numindstack);
varphasered=varphaseall;

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

figure(1);hold on;
param=tau;
scatter(param,varpreind_hanssen,50,'gs');
scatter(param,varprered_hanssen,50,'gd','filled')
scatter(param,varpreind_agram,50,'bs');
scatter(param,varprered_agram,50,'bd','filled')
scatter(param,varpreind_yujie,50,'rs');
scatter(param,varprered_yujie,50,'rd','filled')
end
figure(1)
legend('Hanssen independent','Hanssen redundant', 'Agram independent','Agram redundant','Yujie independent','Yujie redundant')
% xlabel('\rho_\infty')
xlabel('\tau, days')
ylabel('\sigma^2_{predict}')
grid on
fig=gcf;fig.InvertHardcopy = 'off';
saveas(gcf,'varydecor_compare3models','epsc')

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