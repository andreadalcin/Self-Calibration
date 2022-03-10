function [Data,G] = load_h155_kitty( sequencename )
%LOAD_DATA_H155 Summary of this function goes here
%  out=1 load outliers
%  out=0 load inliers only



load(strcat(sequencename,'_truth.mat'));
nPoints = size(y,2);
nFrames = size(y,3);

Data.nAllPoints = nPoints;
Data.yAll = y;
Data.visibleAll = true(nPoints, nFrames);
Data.ySparse = y;
Data.visibleSparse = true(nPoints,nFrames);
Data.SparseIndex = 1:nPoints;
Data.GtLabel = s(:);
Data.nFrames = nFrames;
Data.nSparsePoints = nPoints;
G = s(:);

%%




end

