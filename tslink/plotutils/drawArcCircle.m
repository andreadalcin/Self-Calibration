function [] = drawArcCircle(center,radius, thetaStart, thetaEnd, h,col, lineWidth)

if(nargin <3)
    h = hggroup('Parent', gca);
end
if(nargin < 6)
    col = 'k';
    lineWidth = 0.01;
end
if(nargin < 7)
    lineWidth = 0.01;
end

theta = linspace(thetaStart,thetaEnd,20);
costheta = cos(theta);
sintheta = sin(theta);


x = center(1)+ radius.*costheta;
y = center(2)+ radius.*sintheta;


% Draw the thinner foreground colored circle
line(x,y,'Parent',h, ...
    'Color',col, ...
    'LineWidth', lineWidth);



end

