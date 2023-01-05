function  g=GfunS(u,v,par)
b=par(2); gamma=par(3);
dy2=b-u.^2.*v;
g= gamma*dy2;