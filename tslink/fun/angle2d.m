function angle = angle2d(u,v)
% compute the angle in radiants between 2d vectors
assert(numel(u)==2 && numel(v)==2, "Dimension mismatch, angle should be between 3d vectors");
x1 = u(1);
y1 = u(2);
x2 = v(1);
y2 = v(2);
angle = atan2(x1*y2-y1*x2,x1*x2+y1*y2);
end