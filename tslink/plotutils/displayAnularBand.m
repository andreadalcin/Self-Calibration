function [  ] = displayAnularBand(X, M, epsi , col)
% DISPLAYANULARBAND show a circular segment togheter with its inlier band
assert(epsi>=0, "The inlier threhsold must be positive");
if(nargin==2)
    col = 'r';
end

V = bsxfun(@minus, X, M(1:2));
theta = bsxfun(@atan2, V(2,:),V(1,:));

[theta_min, u_min] = min(theta);
[theta_max, u_max] = max(theta);

if(theta_min<-pi/2 && theta_max>pi/2)
    %     v1 = V(:,u_min);
    %     v2 = V(:,u_max);
    %     angle = angle2d(v1,v2);
    %     if(angle <pi)
    %         theta(theta<0) = theta(theta<0)+2*pi;
    %         theta_min = min(theta);
    %         theta_max = max(theta);
    %
    %     end
    tol = pi/10;
    if(abs(theta_min +2*pi -theta_max) < tol)
        theta(theta<0) = theta(theta<0)+2*pi;
        theta_min = min(theta);
        theta_max = max(theta);
    end
    
end

radius_out = M(3) + epsi;
radius_inn = M(3) - epsi;


k = 100;
t = linspace(theta_min, theta_max,k);
x_inn = radius_out*cos(t) + M(1);
y_inn = radius_out*sin(t) + M(2);
x_out = radius_inn*cos(t) + M(1);
y_out = radius_inn*sin(t) + M(2);
px = [x_inn(1), x_out,flip(x_inn(2:end))];
py = [y_inn(1), y_out,flip(y_inn(2:end))];
patch(px, py, col, 'FaceAlpha',0.2,'EdgeColor',col)
%rectangle('Position',[par(1) - radius_inn, par(2) - radius_inn, radius_inn*2, radius_inn*2],'Curvature',[1,1],'FaceColor','w','EdgeColor',col);
axis square;
end
