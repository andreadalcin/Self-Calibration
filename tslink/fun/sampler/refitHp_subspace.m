function S =  refitHp_subspace(S, X,epsi, opts)
if(~isfield(opts, 'x84'))
    opts.x84 = false;
end
m = size(S.R,2);
dim_subspace = opts.dim_subspace;
cardmss = dim_subspace;
Hnew = S.H;
Rnew = S.R;

maxIter = 2;
for j = 1:m
    for iter = 1:maxIter
        inliers = Rnew(:,j)<epsi;
        if(opts.x84)
            epsix84 = min(epsi, x84_( Rnew(inliers,j), cardmss));
            inliers = Rnew(:,j)<epsix84;
        end
        if(sum(inliers)> cardmss)        
            hRefitted = fit_flat(X(:,inliers),dim_subspace);
            Hnew(:,j) = flatStruct2Vect(hRefitted);
            Rnew(:,j) = res_flat(X,hRefitted);
        end
    end
end
S.H = Hnew;
S.R = Rnew;

end