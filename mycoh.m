function coh=mycoh(slc1,slc2,numlook)
 % return a single number indicating absolute value of the coherence
 igramsl=(slc1.*conj(slc2));
 igram=cpxlooks(igramsl,numlook);
 intensity1=slc1.*conj(slc1);
 intensity1=cpxlooks(intensity1,numlook);
 intensity2=slc2.*conj(slc2);
 intensity2=cpxlooks(intensity2,numlook);
 coh=abs(igram)./sqrt(intensity1)./sqrt(intensity2);
 coh=mean(coh);
 return