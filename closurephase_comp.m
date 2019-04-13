function closurephase=closurephase_comp(numlook,SLCs,i,j,k)
igramij=SLCs(i,:).*conj(SLCs(j,:));igramij=cpxlooks(igramij,numlook);
igramik=SLCs(i,:).*conj(SLCs(k,:));igramik=cpxlooks(igramik,numlook);
igramjk=SLCs(j,:).*conj(SLCs(k,:));igramjk=cpxlooks(igramjk,numlook);

closurephase=igramij.*igramjk.*conj(igramik);
closurephase=angle(closurephase);
figure;
histogram(closurephase(:));
title(['phase mean=' num2str(mean(closurephase(:))) ' std=' num2str(std(closurephase(:)))]);
return