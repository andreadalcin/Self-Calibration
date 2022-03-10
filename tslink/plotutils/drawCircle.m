function [] = drawCircle(center,radius,h, costheta,sintheta, col)

lineWidth = 2;
if(nargin <3)
    h = hggroup('Parent', gca);
end
if(nargin <4)
    thetaResolution = 2;
    theta=(0:thetaResolution:360)'*pi/180;
    costheta = cos(theta);
    sintheta = sin(theta);
end
if(nargin<5)
    col = 'k';
end



x = center(1)+ radius.*costheta;
y = center(2)+ radius.*sintheta;



% Draw the thinner foreground colored circle
line(x,y,'Parent',h, ...
    'Color',col, ...
    'LineWidth', lineWidth);



end

