% simulate closure phase for multilook igrams
close all; clc; clear;

numlooks=10;
%% distributed scatterers
mu=0;sigma=1;
simudata=normrnd(mu,sigma,[3,2000]);
real=simudata(:,1:2:end);imag=simudata(:,2:2:end);
SLCs1=complex(real,imag);
% single look test
igram12=SLCs1(1,:).*conj(SLCs1(2,:));
igram23=SLCs1(2,:).*conj(SLCs1(3,:));
igram13=SLCs1(1,:).*conj(SLCs1(3,:));
triplet=igram12.*igram23.*conj(igram13);
triplet=reshape(triplet,100,10);
phasetriplet=angle(triplet);
figure;imagesc(phasetriplet);colorbar;
title('Single-loook,distributed scatterers')
% 10 looks
igram12_l10=cpxlooks(igram12,numlooks);
igram23_l10=cpxlooks(igram23,numlooks);
igram13_l10=cpxlooks(igram13,numlooks);
phasetriplet=igram12_l10.*igram23_l10.*conj(igram13_l10);
phasetriplet=reshape(phasetriplet,10,10);
figure;imagesc(angle(phasetriplet));colorbar;
title('10-looks,distributed scatterers')
%% point scatterers
mu=0;sigma=1;
simudata=normrnd(mu,sigma,[3,2000]);
simudata(:,1:500)=simudata(:,1:500)+ones(3,500)+rand(3,500);
simudata(:,501:1000)=simudata(:,501:1000)+2.*ones(3,500)+rand(3,500);
simudata(:,1001:2000)=simudata(:,1001:2000)+3.*ones(3,1000)+rand(3,1000);
real=simudata(:,1:2:end);imag=simudata(:,2:2:end);
SLCs2=complex(real,imag);
igram12=SLCs2(1,:).*conj(SLCs2(2,:));
igram23=SLCs2(2,:).*conj(SLCs2(3,:));
igram13=SLCs2(1,:).*conj(SLCs2(3,:));
% 10 looks
igram12_l10=cpxlooks(igram12,numlooks);
igram23_l10=cpxlooks(igram23,numlooks);
igram13_l10=cpxlooks(igram13,numlooks);
triplet=igram12_l10.*igram23_l10.*conj(igram13_l10);
triplet1=triplet(1:25);triplet2=triplet(26:50);triplet3=triplet(51:100);
triplet=reshape(triplet,10,10);
phasetriplet=angle(triplet);
figure;imagesc((phasetriplet));colorbar;
title('10 looks,point scatterers, amp=1,2,3')
phasetriplet1=angle(triplet1);phasetriplet2=angle(triplet2);phasetriplet3=angle(triplet3);
mean(phasetriplet1(:))
std(phasetriplet1(:))
mean(phasetriplet2(:))
std(phasetriplet2(:))
mean(phasetriplet3(:))
std(phasetriplet3(:))
% 
%% mixed
SLCs3=[SLCs1(:,1:500) SLCs2 SLCs1(:,501:1000)];
igram12=SLCs3(1,:).*conj(SLCs3(2,:));
igram23=SLCs3(2,:).*conj(SLCs3(3,:));
igram13=SLCs3(1,:).*conj(SLCs3(3,:));
% 10 looks
igram12_l10=cpxlooks(igram12,numlooks);
igram23_l10=cpxlooks(igram23,numlooks);
igram13_l10=cpxlooks(igram13,numlooks);
phasetriplet=igram12_l10.*igram23_l10.*conj(igram13_l10);
phasetriplet=reshape(phasetriplet,10,20);
figure;imagesc(angle(phasetriplet));colorbar;
title('10 looks,mixed scatterers')