function [ distFun, hpFun, fit_model, cardmss, isdegen, d  ] = set_model( model )
%SET_MODEL set the model specification: possible models are 'line','segment', 'circle',
% fundamental matrices ('fundamental'), 'subspace4', affspace3.
%
% OUTPUT:
% distFun: evaluates the distances between points and model,
% hpFun: genrate the hypothesis function using
% points with minimal sample cardinality (cardmss),
% fit_model: fits models to set of points with cardinality equal or larger than
% cardmss: cardinality of minimal sample set (e.g 2 for lines)

switch model
    case 'line'
        distFun = @distPointLine;
        hpFun = @hpLines;
        fit_model = @fit_line;
        isdegen = @isdummy_degen;
        cardmss = 2;
        d=3;
    case 'segment'
        distFun = @distPointSeg;
        hpFun = @hpSeg;
        fit_model = @fit_lines;
        cardmss = 2;
        d=4;
        
        isdegen= @dummydegenerate;
    case 'circle'
        distFun = @distPointCircle;
        hpFun = @hpCircles;
        fit_model = @fit_circles;
        cardmss = 3;
        isdegen= @isdegen_line;
        d=3;
    case 'fundamental'
        distFun = @distPointFm;
        hpFun = @hpFundamental;
        fit_model = @fit_Fundamental;
        cardmss = 8;
        
        isdegen=@isdummy_degen;%;%
        d=9;
    case 'homography'
        distFun =@distPointH;%@distPointHgeometric;%;%@distPointHgeometric; %%
        hpFun = @hpHomography;
        fit_model = @fit_Homography;
        cardmss = 4;
        d=9;
        isdegen=@isdummy_degen;%@homography_degen;%%@islinedegenerate;% isdegen_line;%@isdummy_degen;%@%@isHomographydegenerate;%@ isdegen_line;%
        %disp('Da fare controllo quadrato')
        
    case 'subspace4'
        distFun = @distPointSubspace4;
        hpFun = @hpSubspace4;
        fit_model = @fit_Subspace4;
        cardmss = 4;
        isdegen= @dummydegenerate;%@isSub4degenerate;
        d = nan;
    case 'affspace3'
        distFun=@distPointAffspace;
        hpFun=@hpAffspace;
        fit_model=@fit_aff;
        cardmss=4;
        isdegen= @dummydegenerate;%@isAff3degenerate;
        d=nan;
    case 'plane'
        distFun = @distPointPlane;
        hpFun = @hpPlanes;
        fit_model = @fit_Planes;
        cardmss = 3;
        d=1;
        isdegen= @dummydegenerate;
        
        
    case 'subspace9'
        distFun = @distPointSubspace9;
        hpFun = @hpSubspace9;
        fit_model = @fit_Subspace9;
        cardmss = 9;
        isdegen= @dummydegenerate;
        d = nan;
        
    case 'trifocal'
        distFun = @distPointTrifocal;
        hpFun = @hpTrifocal;
        cardmss = 7;
        fit_model = @fit_trifocal;
        isdegen =@dummydegenerate;
        d = 27;
    otherwise
        warning('model not yet supported: possible models are line, segment,circle, homography, fundamental, subspace4')
end



end

