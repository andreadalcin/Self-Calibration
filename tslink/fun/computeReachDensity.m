function [lrd] = computeReachDensity(RD,knn)
%COMPUTEREACHDENSITY
n = size(RD,1);
assert(size(RD,2)==n);
assert(numel(knn)==n);
lrd = inf(n,1);
for i =1:n
    avgRD = 0;
    for j = knn{i}
        avgRD = avgRD+ RD(i,j);
    end
    avgRD = avgRD/numel(knn{i});
    lrd(i) = avgRD;
end
lrd = 1./lrd;
end

