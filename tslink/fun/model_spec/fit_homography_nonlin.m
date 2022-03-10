function h = fit_homography_nonlin(X,h_init)
h_init = reshape(h_init,[3,3]);
[h,~] = vgg_H_from_x_nonlin(h_init, X(1:3,:),X(4:6,:));
%h = homog_nonlin(h_init, X(1:2,:), X(4:5,:));
h = h(:);
end