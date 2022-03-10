
function Xp = DataKanatani(X,r)
n=size(X,2);
centroid = mean(X,2);
% deviations
P = bsxfun(@minus, X, centroid);

[U,~,~]= svd(P,0);
U=U(:,1:r);

Xp = U' * P;
end
