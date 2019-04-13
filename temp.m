clear all;clc;close all;
gama=0.3;snr=1;r=1./sqrt(snr);

N=100;
y13=0;y24=0;y14=0;y23=0;
mu=0;sigma=sqrt(0.5);
signal13=0;signal24=0;signal14=0;signal23=0;
noise13=0;
coef1=sqrt(1-gama^2);
for i=1:N
    
    a1=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    d2=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    d3=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    d4=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n1=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n2=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n3=complex(normrnd(mu,sigma),normrnd(mu,sigma));
    n4=complex(normrnd(mu,sigma),normrnd(mu,sigma));

    
    noise=a1*(coef1*(gama*conj(d2)+conj(d3))+r*conj(n3))+r*n1*gama^2*conj(a1);
    noise13=noise13+noise;
    signal13=signal13+gama^2*abs(a1)^2;
    y13=y13+gama^2*abs(a1)^2+noise;
    
    figure(1);hold on;compass(signal13,0,'r');compass(real(noise),imag(noise));title('y13')
    
    a2=gama*a1+coef1*d2;
    noise=a2*(coef1*(gama*conj(d3)+conj(d4))+r*conj(n4))+r*n2*gama^2*conj(a2);
    signal24=signal24+gama^2*abs(a2)^2;
    y24=y24+gama^2*abs(a2)^2+noise;
    figure(2);hold on;compass(signal24,0,'r');compass(real(noise),imag(noise));title('y24')
    
    noise=a1*(coef1*(gama^2*conj(d2)+gama*conj(d3)+conj(d4))+r*conj(n4))+r*n1*gama^3*conj(a1);
    signal14=signal14+gama^3*abs(a1)^2;   
    y14=y14+gama^3*abs(a1)^2+noise;
    figure(3);hold on;compass(signal14,0,'r');compass(real(noise),imag(noise));title('y14')
    
    noise=a2*(coef1*conj(d3)+r*conj(n3))+r*n2*gama*conj(a2);
    signal23=signal23+gama*abs(a2)^2;
    y23=y23+gama*abs(a2)^2+noise;
    figure(4);hold on;compass(signal23,0,'r');compass(real(noise),imag(noise));title('y23')
end
hold off;
figure(5);hold on;
compass(real(y14),imag(y14),'r');
compass(real(y24),imag(y24),'b');
compass(real(y13),imag(y13),'y');
compass(real(y23),imag(y23),'k');
hold off;
legend('y14','y24','y13','y23')
angle_13=angle(y13);
angle_24=angle(y24);
(angle_13+angle_24)/2
 angle_14=angle(y14);
 angle_23=angle(y23);
 (angle_13+angle_24+angle_14+angle_23)/4