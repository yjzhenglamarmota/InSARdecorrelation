function stat_SLC_timecoh(SLCs,numlook)
% scatter plot of time span vs coh
numSLC=size(SLCs,1);
k=1;
numigram=numSLC/2*(numSLC-1);
t=zeros(numigram,1);coh=zeros(numigram,1);
for i=1:numSLC-1
    for j=i+1:numSLC
        t(k)=j-i;
        coh(k)=mycoh(SLCs(i,:),SLCs(j,:),numlook);
        k=k+1;
    end
end
figure;scatter(t,coh,10,'r*')  ;
xlabel('temporal baseline, unitless')
ylabel('coherence')
grid on;

return