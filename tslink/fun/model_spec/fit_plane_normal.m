function [ par ] = fit_plane_normal( X )
%FIT_PLANE using a point  and a normal
% 2021-07
num_points = size(X,2);
if(num_points == 1)
    par = fitplane_point_and_normal(X(1:3,:),X(4:6,:));
elseif(num_points==3)
    par =  fitplane_cardmss(X');
elseif(num_points >3)
    par = fitplane(X);
    par = par./norm(par(1:3));
end
par = par(:);
end

%% fit a plane using a point and a normal
function model = fitplane_point_and_normal(point, normal)
% A(x-x0)+B(y-y0)+C(z-z0)=0
% Ax + By + Cz - (Ax0 + By0 + Cz0)
denom = sum(normal.^2);
if denom < eps(class(point))
    model = [];
else
    normal = normal / sqrt(denom);
    d = -point * normal';
    model = [normal, d];
end
end

%% for cardmss
function model = fitplane_cardmss(points)

a = points(2, :) - points(1, :);
b = points(3, :) - points(1, :);
% Cross product
normal = [a(2).*b(3)-a(3).*b(2), ...
    a(3).*b(1)-a(1).*b(3), ...
    a(1).*b(2)-a(2).*b(1)];
denom = sum(normal.^2);
if denom < eps(class(points))
    model = [];
else
    normal = normal / sqrt(denom);
    d = -points(1,:) * normal';
    model = [normal, d];
end
end
%%
function B = fitplane(XYZ)

[rows,npts] = size(XYZ);

if rows ~=3
    error('data is not 3D');
end

if npts < 3
    error('too few points to fit plane');
end

% Set up constraint equations of the form  AB = 0,
% where B is a column vector of the plane coefficients
% in the form   b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0.

A = [XYZ' ones(npts,1)]; % Build constraint matrix

if npts == 3             % Pad A with zeros
    A = [A; zeros(1,4)];
end

[u d v] = svd(A);        % Singular value decomposition.
B = v(:,4);              % Solution is last column of v.
end