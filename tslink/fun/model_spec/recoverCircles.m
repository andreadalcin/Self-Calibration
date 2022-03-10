function [M,R] = recoverCircles(X,G)
% RECOVERCIRCLES given a labeling of the data G recover the models
% Optionally produce the residuals matrix
assert(size(X,2)==numel(G),'Dimension mismatch between numper of points and labeling');
k = max(G);
M = nan(3,k);
for i = 1:k
    inliers = X(:,G==i);
    if(sum(G==i)>=3)
        par = fit_circle_taubin(inliers);
        % par = fit_circle_lm(inliers,par);
        M(:,i) = par;
    else
        fprintf('there are too few points in cluster %i to define a circle\n',i)
    end
    
end

if(nargout>1)
    n = size(X,2);
    R = nan(n, k);
    for j = 1:k
        R(:,j) = res_circle(X,M(:,j));
    end
end

