function [M,res, thresh] = recover_affine_fundamental_geo( X,G )
n = max(G);
M = nan(9,n);
if(nargout>1)
    res = nan(size(G));
     thresh = nan(1,max(G));
end
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=4)
        M(:,i) = fit_fm_affine(inliers);
        res(G==i) = res_fm_geo(inliers,M(:,i));
        %M(:,i)= fit_fm_bnb(inliers);
        %M(:,i) = fit_fm_andrea(inliers,M(:,i));
        %M(:,i) = fit_fm_torr(inliers,M(:,i));
      
        if(nargout>1)
            u = res_fm(inliers,M(:,i));
            res(G==i) = u;
            thresh(i) = x84_(u,1, 3.5);
        end
    else
        fprintf('there are too few points in cluster %i to define a fm\n',i)
    end

end


end

