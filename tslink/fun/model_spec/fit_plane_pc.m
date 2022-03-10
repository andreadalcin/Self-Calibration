function [m] = fit_plane_pc(X,epsi)
%FIT_PLANE_PC Summary of this function goes here
%   Detailed explanation goes here
 pc = pointCloud(X(1:3,:)');
 model = pcfitplane(pc,epsi);
 m = model.Parameters(:);
end

