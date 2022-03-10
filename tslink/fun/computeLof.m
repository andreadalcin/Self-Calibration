function [lof, lrd, RD, knn] = computeLof(D,k)
%COMPUTELOF compute local outlier factor

n = size(D,1);
assert(n == size(D,2),'Dimension mismatch, distance matrix must be square');
assert(k>0,'k must be a non negative integer');
[RD,knn] = computeReachDist(D,k);
[lrd] = computeReachDensity(RD,knn);

lof = nan(n,1);
for i =1:n
    avgLRD = 0;
    for j = knn{i}
        avgLRD = avgLRD+ lrd(j);
    end
    avgLRD = avgLRD/numel(knn{i});
    lof(i) = avgLRD/lrd(i);
end


end

