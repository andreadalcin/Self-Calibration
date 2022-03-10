function [Xgridded, E] = getAxAlignedGrid(X, row, col)
%getAxAligendGrid reorder the points in a axis aligned rowmajor grid for plotting purposes
% Xgirdded contain the same points of X reordered for plotting purpose.
% E contains the indices of the edges of the grid.
n = size(X,2);
assert(n ==row*col,'Dimension mismatch between number of points and grid layout.');
xx = unique(X(1,:));
yy = unique(X(2,:));
assert(numel(xx)==col,'Points are not in an axed aligned grid.');
assert(numel(yy)==row,'Points are not in an axed aligned grid.');
Xgridded = [];
for y =sort(yy,'descend')
    for x = xx
        Xgridded=[Xgridded, [x;y]];
    end
end

% get edges
E = [];
for k = 1:n
    i =  ceil(k/col);
    j = rem(k-(i-1)*col,col+1);
    if(j < col)
        E = [E; [k, (k +1)]];
    end
    if(i < row)
        E = [E; [k, (k +col)]];
    end
end

end

