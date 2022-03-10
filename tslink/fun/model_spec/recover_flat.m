function [M,res, mads,sigmas] = recover_flat( X,G ,dime)
n = max(G);
M = [];
mads = nan(1,n);
sigmas = nan(1,n);
if(nargout>1)
    res = nan(size(G));
end
for i = 1:n
    inliers = X(:,G==i);
    if(sum(G==i)>=dime+1)
        m = flatStruct2Vect(fit_flat(inliers,dime));
        M = [M,m];
         if(nargout>1)
            u = res_flat(inliers,m);
            res(G==i) = u;
            mads(i) = x84_(u, 1);
            sigmas(i) = std(u);
        end
    else
        fprintf('there are too few points in cluster %i to define a line\n',i)
    end

end


end

