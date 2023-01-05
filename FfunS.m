function  f=FfunS(u,v,par) 
%par=[a, b, gamma];
a=par(1); gamma=par(3);
dy1=a-u+u.^2.*v;
f=gamma*dy1;
