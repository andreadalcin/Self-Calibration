function model = fit_pc_cylinder_ls_circle( X, W )
% project points onto the plane orthogonal to w0;
W = normc(W);
P = eye(3) - W*W';
Y = P*X(1:3,:);
% least square fitting of cylinder;
Y = P*X(1:3,:);
Y = Y([1,3],:);
model_circle = fit_circle_taubin(Y);
a = model_circle(1);
b = model_circle(2);
r = model_circle(3);
PC = [a;0;b];
model = [PC; W; r];
end