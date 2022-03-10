function [M,res] = recover_plane( X,G )
n = max(G);
M = nan(4,n);
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=3)
        M(:,i) = fit_plane(inliers);
    else
        fprintf('there are too few points in cluster %i to define a line',i)
    end
    if(nargout>1)
        res(G==i) = res_plane(inliers,M(:,i));
    end
end


end

