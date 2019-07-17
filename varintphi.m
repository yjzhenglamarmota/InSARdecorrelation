function varphi=varintphi(rho)
% function to compute var(\phi) when L=1
% phi here is interferometric phase
varphi=pi^2/3-pi*asin(rho)+(asin(rho)).^2-0.5*dilog(1-rho.^2);

return
