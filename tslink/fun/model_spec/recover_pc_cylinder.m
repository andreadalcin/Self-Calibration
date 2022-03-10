function [M,res] = recover_pc_cylinder( X,G )
% return finite cylinder model for plot
n = max(G);
M = nan(7,n);
res = nan(numel(G),n);
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=2)
        M(:,i) = fit_pc_cylinder_ls(inliers(1:3,:));
        M(:,i) = convertToFiniteCylinder(M(:,i),X);
    else
        fprintf('there are too few points in cluster %i to define a cylinder',i)
    end
    if(nargout>1)
        res(G==i) = res_pc_cylinder(inliers,M(:,i));
    end
end


end

