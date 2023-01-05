% Snapshots in I1
data.u0 = SnapU(:,1);
data.v0 = SnapV(:,1);
if model == 3
    snap_u = SnapU(:,1:2:index);
    snap_v = SnapV(:,1:2:index);
elseif model == 2
    snap_u = SnapU(:,1:4:index);
    snap_v = SnapV(:,1:4:index);
else
    snap_u = SnapU(:,1:inew);
    snap_v = SnapV(:,1:inew);
end
data.nt = index;
Fn = F(snap_u,snap_v,data.par);
Gn = G(snap_u,snap_v,data.par);

% POD bases in I1
disp('SVD 1')
[data.U,~,~] = svd(snap_u,'econ');
if model == 1
    clear snap_u
end
disp('SVD 2')
[data.V,~,~] = svd(snap_v,'econ');
if model == 1
    clear snap_v
end
disp('SVD 3')
[data.Fnl,~,~] = svd(Fn,'econ');
if model == 1
    clear Fn
end
disp('SVD 4')
[data.Gnl,~,~] = svd(Gn,'econ');
if model == 1
    clear Gn
end
disp('SVD done')

% Set correction values
if model == 2
    R = 35;
    r = 10;
elseif model == 3
    R = 41;
    r = 10;
else 
    r = 10;
    R = 45;
end
% Compute the correction terms
[data.uR,data.vR,time_c] = pod_classic(model,1,R,R,data);
% Compute the solution in I1 by PODc
[u_podc,v_podc,time] = pod_correction(model,1,r,r,R,data);
u_podc = data.U(:,1:r)*u_podc(:,end);
v_podc = data.V(:,1:r)*v_podc(:,end);
time_i1 = time;
