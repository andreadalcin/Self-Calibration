function resi = res_parabola(X,h)
m = size(h,2);
n = size(X,2);
resi = nan(n,m);
for j = 1:m
    H = h(:,j);
    for i = 1:n
        P =X(:,i);
        resi(i,j) = distPointParabola(P,H);
    end
end
end