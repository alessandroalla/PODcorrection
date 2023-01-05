function  [ur,vr,time_podcorr] = pod_correction(model,r,R,data)

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
        u1 = Upod'*data.A*UR*data.uR;
        v1 = Vpod'*data.A*VR*data.vR;
        time_podcorr = 0;
        tic;
        for k = 1 : data.nt-1
            ur(:,k+1) = ur(:,k) + data.ht*(u1(:,k+1)+Upod'*F(UR*data.uR(:,k),VR*data.vR(:,k),data.par));
            vr(:,k+1) = vr(:,k) + data.ht*(data.d.*v1(:,k+1)+Vpod'*G(UR*data.uR(:,k),VR*data.vR(:,k),data.par));
        end
        time_podcorr = toc;

end




     

