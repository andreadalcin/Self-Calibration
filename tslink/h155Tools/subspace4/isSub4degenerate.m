function  y  = isSub4degenerate( X, theta , cardmss )
%ISHDEGENERATE Check if a sample of 4 points X is degenerate i.e <X> <4
%   INPUT:
%          X: fx4 matrix-sample of four points in homogeneous coordinate
%          in two views
%          y: boolean- true if X is H-degenerate.
%          theta: parameter vectors of subspace matrices in columns

%%
tol=1e-2;
m=size(theta,1);
k=m/4;
T=reshape(theta,k,4);
[~,lambda,~]=svd(T);
L=diag(lambda);
y=any(L<tol);
