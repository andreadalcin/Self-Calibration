function [xl,yl] = get2Daxlim(X,k)
%GET2DBB Summary of this function goes here
%   Detailed explanation goes here
if(nargin<2)
    k = 0.01;
end
xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));
wx = xmax-xmin;
wy = ymax -ymin;
dx = k*wx;
dy = k*wy;
xl = [xmin-dx, xmax+dx];
yl = [ymin-dy, ymax+dy];
end

