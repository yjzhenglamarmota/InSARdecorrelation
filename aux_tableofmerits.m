% aux: compute table of merit
clear;clc;close all;
% load DVpredicitons.mat
load CASpredicitons.mat
numpt=length(preind);
erryujie=[preind-varphaseind prered-varphasered];
errpiyush=[preind_piyush-varphaseind  prered_piyush-varphasered];
errhanssen=[preind_hanssen-varphaseind  prered_hanssen-varphasered];

indx1=isfinite(erryujie(:,1)); indx2=isfinite(erryujie(:,2)); 
indx=find(indx1.*indx2);
preind=preind(indx);prered=prered(indx,:);
erryujie=erryujie(indx,:);errpiyush=errpiyush(indx,:);errhanssen=errhanssen(indx,:);
varphaseind=varphaseind(indx); varphasered=varphasered(indx);
varphasedata=mean(sqrt(varphaseind.^2+varphasered.^2));
% centroid, standard deviation 
ct_yujie=sqrt(norm(mean(erryujie,1),2)/varphasedata);
ct_piyush=sqrt(norm(mean(errpiyush,1),2)/varphasedata);
ct_hanssen=sqrt(norm(mean(errhanssen,1),2)/varphasedata);

% % norm 2
% erryujienorm2=sqrt(erryujie(:,1).^2+erryujie(:,2).^2);
% n2_yujie=mean(erryujienorm2);
% 
% errpiyushnorm2=sqrt(errpiyush(:,1).^2+errpiyush(:,2).^2);
% n2_piyush=mean(errpiyushnorm2);
% 
% errhanssennorm2=sqrt(errhanssen(:,1).^2+errhanssen(:,2).^2);
% n2_hanssen=mean(errhanssennorm2);
% 
% % norm 1
% erryujienorm1=abs(erryujie(:,1))+abs(erryujie(:,2));
% n1_yujie=mean(erryujienorm1);
% 
% errpiyushnorm1=abs(errpiyush(:,1))+abs(errpiyush(:,2));
% n1_piyush=mean(errpiyushnorm1);
% 
% errhanssennorm1=abs(errhanssen(:,1))+abs(errhanssen(:,2));
% n1_hanssen=mean(errhanssennorm1);
