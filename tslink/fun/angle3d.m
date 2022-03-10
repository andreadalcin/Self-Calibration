function theta = angle3d(u,v)
%ANGLE2D compute the angle in radians between two 2d vector
assert(numel(u)==3 && numel(v)==3, "Dimension mismatch, angle should be between 3d vectors");
theta = atan2(norm(cross(u,v)),dot(u,v));
end

