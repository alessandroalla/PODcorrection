function [ur,vr,time] = pod_classic_old(model,method,r,rdeim,data)

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

switch method
    case 1 % POD
        tic;
        for k = 1 : data.nt-1
            ur(:,k+1) = Uu\(Lu\(PPu*(ur(:,k) + data.ht*Upod'*F(Upod*ur(:,k),Vpod*vr(:,k),data.par))));
            vr(:,k+1) = Uv\(Lv\(PPv*(vr(:,k) + data.ht*Vpod'*G(Upod*ur(:,k),Vpod*vr(:,k),data.par))));
        end
        time = toc;
    case 2 % DEIM
        Udeim = data.Fnl(:,1:rdeim);
        I = speye(data.n);
        [~,~,indexu] = qr(Udeim','vector');
        indexu = indexu(1:rdeim);
        Vdeim = data.Gnl(:,1:rdeim);
        [~,~,indexv] = qr(Vdeim','vector');
        indexv = indexv(1:rdeim);
        PSIdeim_u = Udeim*pinv(Udeim(indexu,:));
        PSIdeim_v = Vdeim*pinv(Vdeim(indexv,:));
        Cu1 = Upod(indexu,:);
        Cv1 = Vpod(indexu,:);
        Cu2 = Upod(indexv,:);
        Cv2 = Vpod(indexv,:);
        Upod_off = Upod'*PSIdeim_u;
        Vpod_off = Vpod'*PSIdeim_v;
        tic;
        for k = 1 : data.nt-1
            tmpu = ur(:,k)+data.ht*Upod_off*F(Cu1*ur(:,k),Cv1*vr(:,k),data.par);
            tmpv = vr(:,k)+data.ht*Vpod_off*G(Cu2*ur(:,k),Cv2*vr(:,k),data.par);
            ur(:,k+1) = Uu\(Lu\(PPu*tmpu));
            vr(:,k+1) = Uv\(Lv\(PPv*tmpv));
        end
        time = toc;
end
