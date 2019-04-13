gama23=0.2;
% gama13=gama23^2;
% gama24=gama23^2;
% gama14=gama23^3;
gama13=0.1;
gama24=0.1;
gama14=0.08;

var_ind=0.25*(25*exp(-35*gama13)+25*exp(-35*gama24));
var_red=1/16*(25*exp(-35*gama13)+25*exp(-35*gama24)+25*exp(-35*gama14)+25*exp(-35*gama23));