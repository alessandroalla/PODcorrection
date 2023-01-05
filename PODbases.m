disp('SVD 1')
[data.U,~,~] = svd(SnapU,'econ');
disp('SVD 2')
[data.V,~,~] = svd(SnapV,'econ');
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
