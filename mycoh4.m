function coh=mycoh4(slc1,slc2,slc3,slc4,numlook)
 % return a single number indicating absolute value of the coherence
 igramsl=(slc1.*conj(slc2).*conj(slc3).*slc4);
 igram=cpxlooks(igramsl,numlook);
 intensity1=slc1.*conj(slc1);
 intensity1=cpxlooks(intensity1,numlook);
 intensity2=slc2.*conj(slc2);
 intensity2=cpxlooks(intensity2,numlook);
 intensity3=slc3.*conj(slc3);
 intensity3=cpxlooks(intensity3,numlook);
 intensity4=slc4.*conj(slc4);
 intensity4=cpxlooks(intensity4,numlook);
 coh=abs(igram)./sqrt(intensity1)./sqrt(intensity2)./sqrt(intensity3)./sqrt(intensity4);
 coh=mean(coh);
 return