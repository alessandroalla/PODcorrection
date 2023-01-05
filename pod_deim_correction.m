function  [ur,vr,time_podcorr] = pod_deim_correction(model,r,rdeim,R,data)

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
UR = data.U(:,1:R);
VR = data.V(:,1:R);
ur0 = Upod'*data.u0;
vr0 = Vpod'*data.v0;
ur(:,1) = ur0;
vr(:,1) = vr0;

% Solve the corrected reduced model
% uc_t = Upod'*A*UR*uR + Upod'*f(UR*uR,VR*vR)
% vc_t = d*Vpod'*(A*VR*vR) +  Vpod'*g(UR*uR,VR*vR)

        Udeim = data.Fnl(:,1:rdeim);
        I = speye(data.n);
        [~,~,indexu] = qr(Udeim','vector');
        indexu = indexu(1:rdeim);
        Vdeim = data.Gnl(:,1:rdeim);
        [~,~,indexv] = qr(Vdeim','vector');
        indexv = indexv(1:rdeim);
        time_podcorr = 0;
        PSIdeim_u = Udeim*pinv(Udeim(indexu,:));
        PSIdeim_v = Vdeim*pinv(Vdeim(indexv,:));
        cu1 = UR(indexu,:)*data.uR;
        cu2 = UR(indexv,:)*data.uR;
        cv1 = VR(indexu,:)*data.vR;
        cv2 = VR(indexv,:)*data.vR;
        u1 = Upod'*data.A*UR*data.uR;
        v1 = Vpod'*data.A*VR*data.vR;
        Upod_off = Upod'*PSIdeim_u;
        Vpod_off = Vpod'*PSIdeim_v;
        tic;
        for k = 1 : data.nt-1
            ur(:,k+1) = ur(:,k) + data.ht*(u1(:,k+1)+ Upod_off*F(cu1(:,k),cv1(:,k),data.par));
            vr(:,k+1) = vr(:,k) + data.ht*(data.d.*v1(:,k+1)+ Vpod_off*G(cu2(:,k),cv2(:,k),data.par));
        end
        time_podcorr = toc;
end




     

