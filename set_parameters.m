function [par,d,Lx,Nx,Tf,ht,U0,V0] = set_parameters(model)

switch model
    case 1 % FHN
        a = 0.1;
        b = 11;
        gamma = 65.731;
        e = 0.1;
        dc = b*(2*sqrt(1-a)+2-a);
        d = dc*(e^2+1);
        par = [gamma,a,b];
        % Spatial domain
        Lx = pi; Ly = Lx;
        Nx = 101;
        % Initial conditions
        load('RandMatr100.mat');
        U0 = Rs.*1e-03;
        V0 = Rs.*1e-03;
        % Time domain
        Tf = 50;
        ht = 1e-4;
    case 2 % Schnakenberg
        d = 10; 
        a = 0.1; 
        b = 0.9;
        ue = a+b; 
        ve = b/ue^2;
        gamma = 1000;
        par = [a, b, gamma, d, ue, ve];
        % Spatial domain
        Lx = 1; 
        Ly = Lx;
        Nx = 51; 
        % Initial conditions 
        load('RandMatr50.mat');
        U0 = ue + Rs*1e-05;
        V0 = ve + Rs*1e-05;
        % Time domain
        ht = 1e-04; 
        Tf = 2;
    case 3 % DIB
        ro = 25/4; 
        a1 = 10;
        a2 = 1;  
        k2 = 2.5; 
        k3 = 1.5;  
        alpha = 0.5; 
        gamma = 0.2;
        d = 20;
        b = 66; 
        c = 3; 
        dd=c*(1-alpha)*(1-gamma+gamma*alpha)/(alpha*(1+gamma*alpha));
        ue=0; ve=alpha;
        par = [c,b,alpha,gamma,a1,a2,k2,k3,dd,ro,d,ue,ve];
        % Spatial domain
        Lx = 20; Ly = Lx;
        Nx = 101; 
        % Initial conditions
        load('RandMatr100.mat');
        U0 = ue + Rs*1e-05;
        V0 = ve + Rs*1e-05;
        % Time domain
        ht = 1e-3;
        Tf = 100; 
end