function y = f_fhn(u,v,par)
gamma = par(1);
y = gamma*(-u.^3+u-v);
