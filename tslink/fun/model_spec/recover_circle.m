function [M, res, mads, sigmas] = recover_circle( X,G )
n = max(G);
M = nan(3,n);
mads = nan(1,n);
sigmas = nan(1,n);
if(nargout>1)
    res = nan(size(G));
end
for i = 1:n
    
    inliers = X(:,G==i);
    if(sum(G==i)>=3)
        par = fit_circle_taubin(inliers);
        %M(:,i) = fit_circle_lm(inliers,par);
        M(:,i) = par;
         if(nargout>1)
            u = res_circle(inliers,M(:,i));
            res(G==i) = u;
            mads(i) = x84_(u, 1);
            sigmas(i) = std(u);
        end
    else
        fprintf('there are too few points in cluster %i to define a circle\n',i)
    end

end


end

