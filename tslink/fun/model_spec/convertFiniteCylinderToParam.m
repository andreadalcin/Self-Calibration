function [m] = convertFiniteCylinderToParam(modelFinite)
%CONVERTFINITECYLINDERTOPARAM 
% convert a finite cylinder to [w, pc, r]
% w axis direction
% pc point on axis
% r radius
pa = modelFinite(1:3);
pb = modelFinite(4:6);
w = pa-pb;
w = w./norm(w);
pc = pa;
r = modelFinite(7);
m = [w(:); pc(:);r];
end

