function [M,res, mads,sigmas] = recover_fundamental( X,G )
n = max(G);
M = nan(9,n);
mads = nan(1,n);
sigmas = nan(1,n);
if(nargout>1)
    res = nan(size(G));
end
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=8)
        M(:,i) = fit_fm(inliers);
        res(G==i) = res_fm(inliers,M(:,i));
        %M(:,i)= fit_fm_bnb(inliers);
        %M(:,i) = fit_fm_andrea(inliers,M(:,i));
        %M(:,i) = fit_fm_torr(inliers,M(:,i));
      
     if(nargout>1)
            u = res_fm(inliers,M(:,i));
            res(G==i) = u;
            mads(i) = x84_(u, 1);
            sigmas(i) = std(u);
        end
    else
        fprintf('there are too few points in cluster %i to define a fm\n',i)
    end

end


end

