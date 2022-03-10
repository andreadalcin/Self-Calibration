function [M,res, mads,sigmas] = recover_homography( X,G )
n = max(G);
mads = nan(1,n);
sigmas = nan(1,n);
M = nan(9,n);
if(nargout>1)
    res = nan(size(G));
end
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=4)
        M(:,i) = fit_homography(inliers);
        %M(:,i) = fit_homography_nonlin(inliers,M(:,i));
        if(nargout>1)
            u = res_homography(inliers,M(:,i));
            res(G==i) = u;
            mads(i) = x84_(u, 1);
            sigmas(i) = std(u);
        end
    else
        fprintf('there are too few points in cluster %i to define a homography\n',i)
    end

end


end

