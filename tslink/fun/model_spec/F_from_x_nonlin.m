function [F,res] = F_from_x_nonlin(F,m1,m2)

% [F,res] = F_from_x_nonlin(F,m1,m2)
% 
% Compute F using non-linear method which minimizes Sampson's approx to
% geometric reprojection error (called 'Gradient criterion" in Luong,
% Faugeras 1996).
% Returns also the array of distances (non squared).
% An initial estimate of F is required.

% This is a wrapper around P. Torr's function for computing F by minimization
% of the Sampson error with the det(F)=0 constraint.
% The interface is the same as vgg_H_from_x_nonlin (see)


[rm1,cm1]=size(m1);
if (rm1 ~= 3)
    error('This function requires homogeneous coordinates');
end

[rm1,cm1]=size(m2);
if (rm1 ~= 3)
    error('This function requires homogeneous coordinates');
end


no_matches = length(m1);

% force det(F)=0 at input
%[U,S,V] = svd(F);
%S(3,3) = 0;
%F = U*S*V';
%f_init = reshape(F',[9,1]);
f_init = F;
% this is from P. Torr
f = torr_nonlinf_mincon2x2(f_init,m1(1,:),m1(2,:),m2(1,:),m2(2,:), no_matches,1);
F = reshape(f, 3, 3); F = F';

% F_sampson_distance_sqr returns *squared distances
res = sqrt(F_sampson_distance_sqr(F,m1,m2));




%	By Philip Torr 2002
%	copyright Microsoft Corp.
%minimize, subject  to the upper 2x2 norm is 1 as described in Torr and Firzgibbon

% 2011, Fusiello: changed optimset

function f = torr_nonlinf_mincon2x2(f_init, nx1,ny1,nx2,ny2, no_matches, m3)

%make sure it is normalized
f_init = f_init / sqrt(f_init(1)^2 + f_init(2)^2 + f_init(4)^2 + f_init(5)^2);
options = optimset('Display','off','Diagnostics','off', 'Algorithm','active-set');

%the function torr_errf_sse takes as input f and all the extra parameters nx1,ny1,nx2,ny2, m3
f = fmincon(@torr_errf_sse,f_init,[],[],[],[],[],[],@torr_nonlcon_f2x2,options,nx1,ny1,nx2,ny2, m3);
%f = fmincon('torr_errf_sse',f_init,[],[],[],[],[],[],[],options,nx1,ny1,nx2,ny2, m3)


%for nonlinear methods...

function sseC = torr_errf_sse(f,nx1,ny1,nx2,ny2, m3)
%disp('estimating error on f')

%disp('estimating squared errors on f')
f = f /norm(f);

r =    f(1) .* nx1(:).* nx2(:) +   f(2).* ny1(:).* nx2(:) + f(3) .* m3.* nx2(:);
r = r +   f(4) .* nx1(:).* ny2(:) +   f(5) .* ny1(:).* ny2(:)+   f(6) .* m3.* ny2(:);
r = r +   f(7) .* nx1(:).* m3+   f(8) .* ny1(:).* m3+   f(9) .* m3.* m3;
r = r.^2;

fdx1 =  f(1) .* nx2(:) + f(4) .* ny2(:) +   f(7) .* m3;
fdx2 =  f(1) .* nx1(:) + f(2).* ny1(:) + f(3) .* m3;
fdy1 =  f(2).* nx2(:) + f(5) .* ny2(:)+ f(8) .* m3;
fdy2 = f(4) .* nx1(:) + f(5) .* ny1(:)+ f(6) .* m3;

g = (fdx1 .* fdx1 +fdx2 .* fdx2 +fdy1 .* fdy1 +fdy2 .* fdy2);

% for non squared error
%   g = sqrt(fdx1 .* fdx1 +fdx2 .* fdx2 +fdy1 .* fdy1 +fdy2 .* fdy2);
%   g = sqrt(g);

e = r./g;
sseC = norm(e,1); % sum of squares

%	By Philip Torr 2002
%	copyright Microsoft Corp.

function [c,ceq] = torr_nonlcon_f2x2(f, nx1,ny1,nx2,ny2, m3)
%c = ...     % Compute nonlinear inequalities at f.
%ceq = ...   % Compute nonlinear equalities at f.

%g(1) = norm(f) -1.0;
%g(2) = f(1) * (f(5) * f(9) - f(6) * f(8)) - f(2) * (f(4) * f(9) - f(6) * f(7)) + f(3) * (f(4) * f(8) - f(5) * f(7));
c = [];
%what norm should we use!
ceq(1) = sqrt(f(1)^2 + f(2)^2 + f(4)^2 + f(5)^2)- 1;
%ceq(1)= norm(f) -1.0;
ceq(2) = f(1) * (f(5) * f(9) - f(6) * f(8)) - f(2) * (f(4) * f(9) - f(6) * f(7)) + f(3) * (f(4) * f(8) - f(5) * f(7));

