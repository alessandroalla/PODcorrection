function [ur,vr,time] = pod_classic(model,r,data)

switch model
    case 1 %FitzHugh-Nagumo
        F = @f_fhn;
        G = @g_fhn;
    case 2 % Schnakenberg
        F = @FfunS;
        G = @GfunS;
    case 3 % DIB
        F = @FfunD;
        G = @GfunD;
end

Upod = data.U(:,1:r);
Vpod = data.V(:,1:r);
Ar = Upod'*data.A*Upod;
Br = Vpod'*data.A*Vpod;
ur0 = Upod'*data.u0;
vr0 = Vpod'*data.v0;
ur(:,1) = ur0;
vr(:,1) = vr0;
Ir = eye(r);
[Lu,Uu,PPu] = lu(Ir-data.ht*Ar);
[Lv,Uv,PPv] = lu(Ir-data.ht*data.d*Br);


        tic;
        for k = 1 : data.nt-1
            ur(:,k+1) = Uu\(Lu\(PPu*(ur(:,k) + data.ht*Upod'*F(Upod*ur(:,k),Vpod*vr(:,k),data.par))));
            vr(:,k+1) = Uv\(Lv\(PPv*(vr(:,k) + data.ht*Vpod'*G(Upod*ur(:,k),Vpod*vr(:,k),data.par))));
        end
        time = toc;
  
end
