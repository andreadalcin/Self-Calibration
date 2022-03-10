function [ M ] =fit_Subspace4( X, C)
%RANSAC_FIT Summary of this function goes here
%   Detailed explanation goes here
if  ~exist('C','var') || isempty(C)
    [U,~,~]=svd(X);
    U=U(:,1:4);
    M=U(:);
else
    lbl=unique(C);
    N = numel(lbl); %numero di cluster;
    f = size(X,1);
    cardmss = 4;
    M = zeros(f*cardmss,N);
    
    for i=1:N 
        L  = lbl(i);
        points2fit = X(:,C==L);
        [U,~,~]=svd(points2fit);
        U=U(:,1:cardmss);
        M(:,i)=U(:);
    end
end

