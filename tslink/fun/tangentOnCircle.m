function [v, p] = tangentOnCircle(circle,x)
%TANGENTONCIRCLE compute the tangent vecotr of a circle in a given point X
% v is the tangent vector. Optionally return the projection of x on the
% circle.
% Cirlce is a 3 vector collecting center x, y and radius


center = circle(1:2);

% project point onto circle
x = x(:) - center(:); % shift the circle to the origin
theta = atan2(x(2),x(1)); % compute the angle
v = [cos(theta +pi/2); sin(theta+pi/2)];
if(nargout >1)
   p = nan(size(X));
  [p(1),p(2)] = pol2cart(theta,circle(3));
end







end

