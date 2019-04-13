function stat_slcs(SLCs)

intensity=abs(SLCs).^2;phase=angle(SLCs);
% subplot(2,1,1);imagesc(intensity);colorbar;title('intensity');
% subplot(2,1,2);imagesc(phase);colorbar;title('phase');
figure;
subplot(2,1,1);hist(intensity(:));
title(['intensity mean=' num2str(mean(intensity(:))) ' std=' num2str(std(intensity(:)))]);
subplot(2,1,2);hist(phase(:));
title(['phase mean=' num2str(mean(phase(:))) ' std=' num2str(std(phase(:)))]);

return