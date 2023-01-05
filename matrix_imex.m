function [SnapU,SnapV,meanU,incrU,time]= matrix_imex(model,Lx,Nx,Tf,ht,p,par,U0,V0)

% More details about the code and the method
% can be find in 10.1016/j.camwa.2019.10.020

switch model
    case 1 % FitzHugh-Nagumo
        F = @f_fhn;
        G = @g_fhn;
    case 2 % Schnakenberg
        F = @FfunS;
        G = @GfunS;
    case 3 % DIB
        F = @FfunD;
        G = @GfunD;
end

L = Lx;
hs = L/(Nx);
N = Nx;
tol = 1e-08;
d = par(end-2);
T = ECDF_McD(p, Nx-1, hs);
I = eye(Nx-1);
Nt = Tf/ht;
tt = linspace(0,Tf,Nt+1);
A1 = I-ht*T;
B = -ht*T';
A2 = I-d*ht*T;
tic;
[Qa1,Ra1] = eig(A1);  %A=Qa*Ra*inv(Qa)
[Qa2,Ra2] = eig(A2);  %A=Qa*Ra*inv(Qa)
[Qb,Rb] = eig(B);  %B=Qb*Rb*inv(Qb)
I_Qa1 = inv(Qa1);
I_Qa2 = inv(Qa2);
I_Qb = inv(Qb);
dA1 = diag(Ra1);
dA2 = diag(Ra2);
dB = diag(Rb);
Inveig_u = 1./(dA1*ones(1,N-1) + ones(N-1,1)*dB');
Inveig_v = 1./(dA2*ones(1,N-1) + ones(N-1,1)*d*dB');
V2 = V0;
U2 = U0;
U21 = U0;
SnapU(:,1) = U0(:);
SnapV(:,1) = V0(:);
meanU(1) = mean(mean(U0));
count = 1;
for k = 1 : Nt
    C_u = U2 + ht*F(U2, V2, par);
    C_v = V2 + ht*G(U2, V2, par);
    Cc_u = I_Qa1*C_u*Qb;
    Cc_v = I_Qa2*C_v*Qb;
    Xx = Cc_u.*Inveig_u;
    Yy = Cc_v.*Inveig_v;
    U2 = Qa1*Xx*I_Qb;
    V2 = Qa2*Yy*I_Qb;
    meanU(k+1)= mean(mean(U2));
    incrU(k) = norm(U21(:)-U2(:));
    if model == 1
        if rem(k,4) == 0
            count = count +1;
            SnapU(:,count) = U2(:);
            SnapV(:,count) = V2(:);
        end
    else
        SnapU(:,k+1)= U2(:);
        SnapV(:,k+1)= V2(:);
    end
    U21 = U2;
end
time = toc;