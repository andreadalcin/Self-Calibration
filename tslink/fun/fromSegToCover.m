function [cov] = fromSegToCover(seg,R,epsi)
%FROMSEGTOCOVER convert a segmentation to a cover
% the segmentaiton is assumed to have the minum numbers of lables
% outliers are marked as 0.

nClusters = max(seg);
nPoints = numel(seg);
cov = zeros(nPoints, nClusters);

if(nargin == 1)
    
    for j = 1:max(seg)
        inliers = (seg == j);
        cov(inliers, j) = 1;
    end
end

if(nargin >1)
    cov = double(R<epsi);
end

