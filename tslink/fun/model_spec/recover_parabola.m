function M = recover_parabola(X,C)
m = max(C);
M = nan(3,m);
for i = 1:m
    inl = C==i;
    if(sum(inl)>2)
        M(:,i) = polyfit(X(1,inl),X(2,inl),2);
    end
end
end