function [] = displayLineInCircle(line,center, radius,h)
%DISPLAYPREFLINES display a line passing through a circle defined by a center and a
%radius
if(nargin < 4)
   h = hggroup('Parent', gca);
end
a = line(1);
b = line(2);
c = line(3);
tol = 1e-3;
if(abs(b)<tol)
    % vertical line
    slope = Inf;
    intercpt = -c/a;
else
    slope = -a/b;
    intercpt = -c/b;
end

[xout,yout] = linecirc(slope,intercpt,center(1),center(2),radius);
line(xout,yout,'Color','k','Linewidth',0.5,'Parent',h)
end

