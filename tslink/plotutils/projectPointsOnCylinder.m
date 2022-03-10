function [pointsProj]=projectPointsOnCylinder(points,model)
N = model(1:3);
pc = model(4:6);
r = model(7);
% (2) determine an othonormal basis q1,q2:
% q1 is the projection of e3 on the plane
e3 = [0;1;0];
q1 = e3 - (e3'*N)/(N'*N) * N;
q1 = normc(q1);
% q2 is the vector on the plane orthogonal to the normal and q1
q2 = cross(N, q1); 
q2 = q2(:);
q2 = normc(q2);
% (3) compute the coords of points w.r.t q1 q2
pointsProj =nan( size(points));
n = size(points,2);
x = nan(1,n);
y = nan(1,n);
u = nan(3,n);
v = nan(3,n);
for i =1:n
    x(i)=(points(:,i)-pc)'*q1;
    y(i)=(points(:,i)-pc)'*q2;
    u(:,i) = x(i)*q1 + y(i)*q2;
    v(:,i) = points(:,i)-(pc + u(:,i));
    u(:,i) = u(:,i)./norm(u(:,i));
    pointsProj(:,i) = pc + r.*u(:,i) + v(:,i);
    
end
