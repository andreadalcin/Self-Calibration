function [M,res] = recover_subspace( X,G,delta )
n = max(G);
M = cell(1,n);
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=delta)
        M{i} = fit_subspace(inliers,delta);
    else
        disp('there are too few points in cluster %i to define a line',i)
    end
    if(nargout>1)
        res(G==i) = res_subspace(inliers,M{i});
    end
end


end

