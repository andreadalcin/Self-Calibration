function [T] = fit_subspace(X,delta)
%FIT_SUBSPACE fit a subspace of a given dimension delta

[k,n]=size(X);

if(delta >= k)
    error('Dimension requested for the plane equal or greater than the dimension of the data')
end

if(n < delta)
    error('Number of points less than the dimension of the hyperplane')
end

if(n == delta)
%The model is the matrix that gives the component of the point orthogonal
%to the plane
%The configuration of the points in x is assumed to be non-degenerate
 y = gramsmithorth(X);
 P = eye(size(X,1))-y*y';
 T.P = P;
elseif(n > delta)
    [U,~,~]=svd(X);
    B = U(:,1:delta);
    T.P =  eye(size(X,1))-B*B';
end
end

