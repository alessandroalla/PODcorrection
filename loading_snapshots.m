switch model
    case 1 % FitzHugh-Nagumo
        F = @f_fhn;
        G = @g_fhn;
        [data.par,data.d,Lx,Nx,Tf,data.ht,U0,V0] = set_parameters(model);
        nx = Nx-1;
        hs = Lx/Nx;
        data.n = nx^2;
        % Time interval
        t = [0: data.ht : Tf];
        data.nt = length(t);
        [SnapU,SnapV,meanU,incrU,time]= matrix_imex(model,Lx,Nx,Tf,data.ht,2,data.par,U0,V0);
        [val,index] = max([zeros(1,100000) incrU(100001:end)]);
    case 2 % Schnakenberg
        F = @FfunS;
        G = @GfunS; 
        [data.par,data.d,Lx,Nx,Tf,data.ht,U0,V0] = set_parameters(model);
        hs = Lx/Nx;
        nx = Nx-1;
        t = [0: data.ht : Tf];
        [SnapU,SnapV,meanU,incrU,time]= matrix_imex(model,Lx,Nx,Tf,data.ht,2,data.par,U0,V0);
        data.n = nx^2;
        data.nt = size(SnapU,2);
        [val,index] = max(incrU);
    case 3 % DIB
        F = @FfunD;
        G = @GfunD;
        [data.par,data.d,Lx,Nx,Tf,data.ht,U0,V0] = set_parameters(model);
        hs = Lx/Nx;
        nx = Nx-1;
        t = [0: data.ht : Tf];
        [SnapU,SnapV,meanU,incrU,time]= matrix_imex(model,Lx,Nx,Tf,data.ht,2,data.par,U0,V0);
        data.n = nx^2;
        data.nt = size(SnapU,2);
        [val,index] = max(incrU);
end
T = diag(-2*ones(nx,1)) + diag(ones(nx-1,1),1) + diag(ones(nx-1,1),-1);
B = zeros(nx,nx);
B(1,1:2) = [2 -1/2];
B(end,end-1:end) = [-1/2 2];
B = 2/3.*B;
T = T + B;
I = speye(nx);
data.A = 1/(hs^2)*kron(I,sparse(T)) + 1/(hs^2)*kron(sparse(T),I);
