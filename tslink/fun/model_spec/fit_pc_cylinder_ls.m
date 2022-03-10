function model = fit_pc_cylinder_ls( X, w0 )
if(nargin<2)
    w0 = [];
end
% least square fitting of cylinder;
[ W, PC, r ,e] = fit_cylinder( X(1:3,:),w0 );
% C point on axis
% r radius of cylinder
% W axis unit direction

model = [PC; W; r];
end