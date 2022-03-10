function [RD,knn] = computeReachDist(D,k)
%COMPUTEREACHDIST compute the reachability distance
% Takes as input a distance matrix D and the k-nn value k
% returns the rachability distance matrix:
% RD(i,j) = max(knn-dist, D(i,j))
% and for each point the list of its knn

n = size(D,2);
assert(n==size(D,1),'Distance matrix must be a square matrix');
knn = cell(1,n);
RD = Inf(n);

[d_sorted,inds] = sort(D,2); % sort the distances
dknn = d_sorted(:,k);        % get the knn-dist
% get the indices of the nearest neighbour
for i = 1:n
    knn{i} = inds(i,d_sorted(i,:)<=dknn(i));
end
% compute the reachability distance matrix
for i = 1:n
    for j = 1:n
        RD(i,j) = max(dknn(j),D(i,j));
    end
end

end

