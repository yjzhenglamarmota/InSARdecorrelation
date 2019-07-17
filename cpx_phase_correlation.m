% test the correlation between complex numbers and the correlation between their corresponding phases.
clear; close all; clc;

numpixel=40000;
cc=1:-0.01:0 ; numSLC=length(cc);
D=simulateSLC(numSLC,numpixel,cc);
numlook=1;
D_ml=cpxlooks(D,numlook);
numpixel=numpixel/numlook;
phase=angle(D_ml);
phase1=phase(1,:);

ccnmr=zeros(numSLC,1);ccphnmr=zeros(numSLC,1);
varphiuw=zeros(numSLC,1);varphiw=zeros(numSLC,1);
for i=1:numSLC
    phasei=phase(i,:);
    phasedf=phasei-phase1;
    varphiuw(i)=var(phasedf);
    varphiw(i)=var(angle(exp(1j*(phasei-phase1))));
    ccnmr(i)=abs(mean(D(1,:).*conj(D(i,:))));
    ccphnmr1=(mean((phase1-mean(phase1)).*(phasedf-mean(phasedf))));
    ccphnmr(i)=1+ccphnmr1/sqrt(var(phase1).*var(phasei));
end
sz=140;
figure(1);scatter(cc,ccphnmr,sz, 'd',... 
    'MarkerEdgeColor',[0, 0.5, 0.5],'MarkerFaceColor',[0,0.7,0.7],'linewidth',2);grid on;
hold on;
yy=0:0.1:1;xx=yy;
plot(xx,yy,'r--','linewidth',2);
xlim([0,1]);ylim([0,1])

% ok let's try theoretical curve
varphi=varintphi(cc);
varphiuwtheo=2/3*pi^2.*(1-cc).^0.45;
ccphtheo=1-(1-cc).^0.45;
scatter(cc, ccphtheo,sz, 'p','MarkerEdgeColor',[0.5, 0, 0.5],'MarkerFaceColor',[0.7,0,0.7],'linewidth',2)


figure;hold on;
plot(cc,varphiuw); plot(cc,varphiw); plot(cc,varphi);plot(cc,varphiuwtheo)
legend('unwrapped','wrapped','wrapped theory','unwrapped theory')