function S =  refitHp(S, X,epsi, opts)
if(~isfield(opts, 'x84'))
    opts.x84 = false;
end
m = size(S.R,2);
model = opts.model;
if(strcmp(model,'line'))
    fitModel = @fit_line;
    resModel = @res_line;
    cardmss = 2;
elseif(strcmp(model,'circle'))
    fitModel = @fit_circle;
    resModel = @res_circle;
    cardmss = 3;
elseif(strcmp(model,'parabola'))
    fitModel = @fit_parabola;
    resModel = @res_parabola;
    cardmss = 3;
elseif(strcmp(model,'homography'))
    fitModel = @fit_homography;
    resModel = @res_homography_geo;
    cardmss = 4;
elseif(strcmp(model,'affine_fundamental'))
    fitModel = @fit_fm_affine;
    resModel = @res_fm_geo;
    cardmss = 4;
 elseif(strcmp(model,'fundamental'))
    fitModel = @fit_fm;
    resModel = @res_fm_geo;
    cardmss = 8;
elseif(strcmp(model,'plane'))
    fitModel = @fit_plane;
    resModel = @res_plane;
    cardmss = 3;
elseif(strcmp(model,'cylinder'))
    fitModel = @fit_cylinder_pc;
    resModel = @res_pc_cylinder;
    cardmss = 100;
else
    error('unsupported model');
end

Hnew = S.H;
Rnew = S.R;

maxIter = 2;
for j = 1:m
    for iter = 1:maxIter
        inliers = Rnew(:,j)<epsi;
        if(opts.x84)
           epsix84 = x84_( Rnew(inliers,j), cardmss);
           inliers = Rnew(:,j)<epsix84;
        end
        if(sum(inliers)> cardmss)
            if(strcmp(model,'cylinder'))
                 hRefitted = fitModel(X(:,inliers),epsi);
            elseif(strcmp(model,'plane'))
                 hRefitted = fitModel(X(1:3,inliers));
            else
                 hRefitted = fitModel(X(:,inliers));
            end
           
            Hnew(:,j) = hRefitted(:);
            Rnew(:,j) = resModel(X,hRefitted);
        end
    end
end
S.H = Hnew;
S.R = Rnew;

end