function [ d ] = distPointSubspace9( X, U )
%DISTPOINTSUBSPACE4 Summary of this function goes here
%   X punto = traiettoria
%   U sottospazio ortonormale di dimensione 4 


f = size(X,1);

U = reshape(U,f,9);
% proietto  X su U ottenendo Y
k=X'*U; %coefficienti di Y in U
Y = k(1).*U(:,1)+k(2).*U(:,2)+k(3).*U(:,3)+k(4).*U(:,4)+k(5).*U(:,5)+k(6).*U(:,6)+k(7).*U(:,7)+k(8).*U(:,8)+k(9).*U(:,9);
d= norm(X-Y);


