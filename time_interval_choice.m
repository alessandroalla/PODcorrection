switch choice
    case 0 % Time interval [0,T]
        data.u0 = SnapU(:,1);
        data.v0 = SnapV(:,1);
        if model >= 2
            SnapU = SnapU(:,1:4:end);
            SnapV = SnapV(:,1:4:end);
        end
    case 1 % Zone I1 reactivity
        data.u0 = SnapU(:,1);
        data.v0 = SnapV(:,1);
        if model == 2
            SnapU = SnapU(:,1:4:index);
            SnapV = SnapV(:,1:4:index);
        elseif model == 3
            SnapU = SnapU(:,1:2:index);
            SnapV = SnapV(:,1:2:index);
        else
            inew = fix(index/4)+1;
            SnapU = SnapU(:,1:inew);
            SnapV = SnapV(:,1:inew);
        end
        data.nt = index;
    case 2 % Zone I2
        % Compute the solution in I1
        find_sol_i1 
        % Set the solution at the final time tau in I1 as initial condition
        % in I2
        data.u0 = u_podc;
        data.v0 = v_podc;
        % Set snapshots in I2
        if model >=2
            SnapU = SnapU(:,index:4:end);
            SnapV = SnapV(:,index:4:end);
        else
            inew = fix(index/4)+1;
            SnapU = SnapU(:,inew:end);
            SnapV = SnapV(:,inew:end);
        end
    t2 = [t(index):data.ht:Tf];
    data.nt = length(t2);           
end

% Solution at the final time in the chosen time interval
u_tf = SnapU(:,end);
v_tf = SnapV(:,end);
% Nonlinearities snapshots
Fn = F(SnapU,SnapV,data.par);
Gn = G(SnapU,SnapV,data.par);