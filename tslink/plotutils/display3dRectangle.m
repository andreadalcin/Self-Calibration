function [ p, projected_points ] = display3dRectangle( theta, points, colspec )
%DISPLAY_3D_RECTANGLE in order to visualize a plane display a 3d rectangular patches 
% INPUT: theta plane parameters vector
%        points 3d points cloud arranged in a 3xn matrix
% OUTPUT: p planar patches


%***********************
alpha_level = 0.5;
%***********************

if(nargin<3)
    colspec= 'y';
end

%% (1) select an origin on points clouds
Origin = points(:,1);
Origin = Origin(:);
e3 = [0;0;1]; % rectangle are aligned with the y-axis
N = theta(1:3); N=N(:); % normal of the plane


%% (2) determine an othonormal basis q1,q2:
% q1 is the projection of e3 on the plane
q1 = e3 - (e3'*N)/(N'*N) * N;
q1 = normc(q1);

% q2 is the vector on the plane orthogonal to the normal and q1
q2 = cross(theta(1:3), q1); 
q2 = q2(:);
q2 = normc(q2);

%% (3) compute the coords of points w.r.t q1 q2
n = size(points,2);
x = nan(1,n);
y = nan(1,n);

for i =1:n
    x(i)=(points(:,i)-Origin)'*q1;
    y(i)=(points(:,i)-Origin)'*q2;
end

%% (4) calculate the extrema on the q1-q2 plane:
a = min(x);
b = min(y);
c = max(x);
d = max(y);

%  (a,d)---------(c,d)
%    |             |
%    |             |
%    |             |
%  (a,b)---------(c,b)

%% (5) and express the vertices of the rectangle in the word reference frame
v1 = Origin + a*q1 + b*q2;
v2 = Origin + a*q1 + d*q2;
v3 = Origin + c*q1 + d*q2;
v4 = Origin + c*q1 + b*q2;
%% (5) bis compute the projected points
projected_points = nan(size(points));
for i = 1:n
    projected_points(:,i) = Origin+x(i)*q1 +y(i)*q2;
end
%%
%% (6) display

xx=[v1(1), v2(1), v3(1), v4(1)];
yy=[v1(2), v2(2), v3(2), v4(2)];
zz=[v1(3), v2(3), v3(3), v4(3)];


p=patch(xx,yy,zz,colspec);
p
set(p,'FaceAlpha',alpha_level,'EdgeAlpha',1,'EdgeColor',colspec,'LineWidth',2);




end

