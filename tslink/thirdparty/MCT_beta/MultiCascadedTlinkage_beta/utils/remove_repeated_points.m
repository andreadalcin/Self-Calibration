function [ Y,flg ] = remove_repeated_points(X, score)
%REMOVE_REPEATED_POINTS remove data points removing repeated point
%   Detailed explanation goes here
N=size(X,2);
flg=true(1,N);
D=pdist(X','euclidean');
D=squareform(D);
Dx = squareform(pdist(X(1:2,:)','euclidean'));
Dy = squareform(pdist(X(4:5,:)','euclidean'));
tol = 1e-12;

for i=1:N
    for j=i+1:N
        if(D(i,j)<tol)% || Dx(i,j)<tol || Dy(i,j)<tol)
            %if(score(i)<score(j))
            %    flg(i)=false;
            %else
                flg(j)=false;
            %end
            
        end
    end
end
Y=X(:,flg);



end

