clear;clc;close all;
numsimu=1000;
angle_13=zeros(numsimu,1);angle_24=zeros(numsimu,1);
angle_23=zeros(numsimu,1);angle_14=zeros(numsimu,1);
ind=zeros(numsimu,1);redu=zeros(numsimu,1);

for simu=1:numsimu
gama=0.4;snr=10;r=1./sqrt(snr);
N=100;
y13=0;y24=0;y14=0;y23=0;
mu=0;sigma=sqrt(0.5); 
signal13=0;signal24=0;signal14=0;signal23=0;

coef1=sqrt(1-gama^2);
atmphase=4/3*pi*rand(4,1);


for i=1:N
    
    d2=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    d3=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    d4=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n1=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n2=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n3=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n4=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    a1=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    a2=gama*a1+coef1*d2;a3=gama*a2+coef1*d3;a4=gama*a3+coef1*d4;
    y1=a1+r*n1;y2=a2+r*n2;y3=a3+r*n3;y4=a4+r*n4;
   
    y1=y1*exp(1j*atmphase(1));y2=y1*exp(1j*atmphase(2));
    y3=y3*exp(1j*atmphase(3));y4=y4*exp(1j*atmphase(4));
    
    signal13=signal13+gama^2*abs(a1)^2*exp(1j*(atmphase(1)-atmphase(3)));
    y13=y13+y1*conj(y3);noise=y1*conj(y3)-gama^2*abs(a1)^2*exp(1j*(atmphase(1)-atmphase(3)));
    figure(1);
    subplot(2,2,1);hold on;compass(real(signal13),imag(signal13),'r');compass(real(noise),imag(noise));title('y13')
  
    signal24=signal24+gama^2*abs(a2)^2*exp(1j*(atmphase(2)-atmphase(4)));
    y24=y24+y2*conj(y4);noise=y2*conj(y4)-gama^2*abs(a2)^2*exp(1j*(atmphase(2)-atmphase(4)));
    subplot(2,2,2);hold on;compass(real(signal24),imag(signal24),'r');compass(real(noise),imag(noise));title('y24')

    signal14=signal14+gama^3*abs(a1)^2*exp(1j*(atmphase(1)-atmphase(4)));
    y14=y14+y1*conj(y4);noise=y1*conj(y4)-gama^3*abs(a1)^2*exp(1j*(atmphase(1)-atmphase(4)));
    subplot(2,2,3);hold on;compass(real(signal14),imag(signal14),'r');compass(real(noise),imag(noise));title('y14')
 
    signal23=signal23+gama*abs(a2)^2*exp(1j*(atmphase(2)-atmphase(3)));
    y23=y23+y2*conj(y3);noise=y2*conj(y3)-gama*abs(a2)^2*exp(1j*(atmphase(2)-atmphase(3)));
    subplot(2,2,4);hold on;compass(real(signal23),imag(signal23),'r');compass(real(noise),imag(noise));title('y23')   
 end

noise13=y13-signal13;
noise24=y24-signal24;
noise14=y14-signal14;
noise23=y23-signal23;

figure(2);hold on;title('independent stack')
compass(real(y13+y24),imag(y13+y24),'-.b');
compass(real(signal13+signal24),imag(signal13+signal24),'r');
compass(real(noise13+noise24),imag(noise13+noise24),'k')
legend('total','signal','noise')
figure(3);hold on;title('redundant stack')
compass(real(y13+y24+y23+y14),imag(y13+y24+y23+y14),'-.b');
compass(real(signal13+signal24+signal23+signal14),imag(signal13+signal24+signal23+signal14),'r');
compass(real(noise13+noise24+noise23+noise14),imag(noise13+noise24+noise23+noise14),'k')
legend('total','signal','noise')

% figure;
% subplot(2,2,1);hold on;compass(signal13,0,'r');compass(real(noise13),imag(noise13),'k');title('y13');
% compass(real(y13),imag(y13),'-.b');
% subplot(2,2,2);hold on;compass(signal24,0,'r');compass(real(noise24),imag(noise24),'k');title('y24');
% compass(real(y24),imag(y24),'-.b');
% subplot(2,2,3);hold on;compass(signal14,0,'r');compass(real(noise14),imag(noise14),'k');title('y14');
% compass(real(y14),imag(y14),'-.b');
% subplot(2,2,4);hold on;compass(signal23,0,'r');compass(real(noise23),imag(noise23),'k');title('y23');
% compass(real(y23),imag(y23),'-.b');
% saveas(gcf,'look_100_2','png')

y_13(simu)=y13;y_24(simu)=y24;y_14(simu)=y14;y_23(simu)=y23;
angle_13(simu)=angle(y13);angle_24(simu)=angle(y24);
ind(simu)=(angle_13(simu)+angle_24(simu))/2;
% ind(simu)=angle(y13+y24) ;
angle_14(simu)=angle(y14);angle_23(simu)=angle(y23);
redu(simu)=(angle_13(simu)+angle_24(simu)+angle_14(simu)+angle_23(simu))/4;
% redu(simu)=angle(y13+y24+y14+y23);
end
std(ind)
std(redu)

figure;
subplot(2,2,1);hold on;histogram(angle_13);title(['y13 std=' num2str(std(angle_13))]);ylim([0,160]);
subplot(2,2,2);hold on;histogram(angle_24);title(['y24 std=' num2str(std(angle_24))]);ylim([0,160]);
subplot(2,2,3);hold on;histogram(angle_14);title(['y14 std=' num2str(std(angle_14))]);ylim([0,160]);
subplot(2,2,4);hold on;histogram(angle_23);title(['y23 std=' num2str(std(angle_23))]);ylim([0,160]);
saveas(gcf,'look_100_phasenoise','png')

phasenoise=cell(4,1);
phasenoise{1}=angle_13';phasenoise{2}=angle_24';
phasenoise{3}=angle_14';phasenoise{4}=angle_23';

for i=1:4
    for j=i:4
        phnoisecov(i,j)=mycoh(phasenoise{i},phasenoise{j},1000);
        phnoisecov(j,i)=phnoisecov(i,j);
    end
end

mycoh(phasenoise{1}+phasenoise{2},phasenoise{3}+phasenoise{4},1000)

expectedstd_red=[0.25*std(angle_13),0.25*std(angle_24),0.25*std(angle_14),0.25*std(angle_23)]...
    *phnoisecov*[0.25*std(angle_13),0.25*std(angle_24),0.25*std(angle_14),0.25*std(angle_23)]';
expectedstd_ind=[0.5*std(angle_13),0.5*std(angle_24)]*phnoisecov(1:2,1:2)*[0.5*std(angle_13),0.5*std(angle_24)]';
expectedstd_red=sqrt(expectedstd_red);expectedstd_ind=sqrt(expectedstd_ind);

% compute the covariance matrice between the igram complex number.
phasenoise=cell(4,1);
phasenoise{1}=y_13;phasenoise{2}=y_24;
phasenoise{3}=y_14;phasenoise{4}=y_23;
for i=1:4
    for j=i:4
        igramnoisecov(i,j)=mycoh(phasenoise{i},phasenoise{j},1000);
        igramnoisecov(j,i)=igramnoisecov(i,j);
    end
end
mycoh(phasenoise{1}+phasenoise{2},phasenoise{3}+phasenoise{4},1000)
% expectedstd_redcp=[0.25,0.25,0.25,0.25]*igramnoisecov*[0.25,0.25,0.25,0.25]';
% expectedstd_indcp=[0.5,0.5]*igramnoisecov(1:2,1:2)*[0.5,0.5]';
% expectedstd_redcp=sqrt(expectedstd_redcp);expectedstd_indcp=sqrt(expectedstd_indcp);
% figure
% scatter(1:simu,ind-redu)
% title('Difference between two stacks')
% 
figure
subplot(1,2,1);histogram(ind);title('Independent stack phase estimation');xlim([-4,4]);
subplot(1,2,2);histogram(redu);title('Redundant stack phase estimation');xlim([-4,4]);
saveas(gcf,'look_100_ind_redu_histogram','png')
% 
