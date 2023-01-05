function  g=GfunD(eta,theta,par)
C=par(1); gam=par(4);
k2=par(7); k3=par(8); D=par(9); ro=par(end-3);
dy2=C*(1+k2.*eta).*(1-theta).*(1-gam.*(1-theta))-D.*theta.*(1+gam.*theta)- ...
D*k3*eta.*theta.*(1+gam*theta);
g= ro*dy2;