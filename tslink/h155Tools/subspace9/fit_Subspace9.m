function [ M ] =fit_Subspace9( X, C)
%RANSAC_FIT Summary of this function goes here
%   Detailed explanation goes here
if  ~exist('C','var') || isempty(C)
    [U,~,~]=svd(X);
    U=U(:,1:4);
    M=U(:);
else
    label=unique(C);
    N=length(label); %numero di cluster;
    f=size(X,1);
    cardmss=9;
    M=zeros(f*9,N);
    
    for i=1:N
        
        L  = label(i);
        points2fit = X(:,C==L);
        [U,D,V]=svd(points2fit);
        U=U(:,1:9);
        M(:,i)=U(:);
    end
    
end
