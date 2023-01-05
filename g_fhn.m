function y = g_fhn(u,v,par)
gamma = par(1);
a = par(2);
b = par(3);

y = gamma*b*(u-a*v);
