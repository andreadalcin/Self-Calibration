function [ H ] = hpSubspace9( X,S)
%HPSUBSPACE4 generate subspace of dimension 4;
%   IDEA genero le m ipotesi, 
m=size(S,1);
f = size(X,1); % lunghezza della traiettoria: delle colonne di X

H = zeros(9*f,m); 
for i=1:m
    [U,~,~]=svd(X(:,S(i,:)));
    O=U(:,1:9);
    %O = orth(X(:,S(i,:)));
   H(:,i)=O(:); 
end

