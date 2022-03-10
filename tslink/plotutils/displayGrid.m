function [] = displayGrid(X, E ,clr)
%DISPLAYGRID display of a grid (aimed at testing synthetic homographies)
% X is assumed to be ordered with getAxAlignedGrid
lineWidthValue = 3;
if(nargin <3)
    clr = 'k';
end
% plot segments
for e = 1:size(E,1)
    x = [X(1, E(e,1)), X(1,E(e,2))];
    y = [X(2, E(e,1)), X(2,E(e,2))];
    line(x,y,'Color',clr,'LineWidth',lineWidthValue);
end


end

