function  f=FfunD(eta,theta,par) 
%par=[C,B,alpha,gam,A1,A2,k2,k3,D,d,ro];
B=par(2); alpha=par(3); A1=par(5); A2=par(6); ro=par(end-3);
%keyboard
dy1=A1*(1-theta).*eta-...
    A2*eta.*eta.*eta-...%^3-...
    B*(theta-alpha);
f=ro*dy1;