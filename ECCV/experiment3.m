clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

load('Dataset/Fs.mat')

width = 4000;
height = 3000;
f = 2 * (width + height);
u0 = width / 2;
v0 = height / 2;

A = [f 0 u0; 0 f v0; 0 0 1];

%%
num_cameras = 10;
% Ps = cell(nchoosek(10,2),1);
Ps = [];
Vs = cell(nchoosek(10,2),1);
index = 1;

for i = 1:num_cameras-1
    for j = i+1:num_cameras
        % P = vgg_P_from_F(Fs(:,:,i,j,1));
        R = rotYPR(rand(1) * 45, rand(1) * 45, rand(1) * 45);
        P = A * [eye(3), 1e1 * rand(3,1)];
        % visualize(P);
        Ps(:,:,index) = P;
        % Vs{index} = 0.5 * [sqrt(width^2 + height^2), 0, width; 0, sqrt(width^2 + height^2), height; 0, 0, 2];
        index = index + 1;
    end
end

centers = [v0 u0];

Vs = 0.5 * [sqrt(width^2 + height^2), 0, width; 0, sqrt(width^2 + height^2), height; 0, 0, 2];
[rPPMs H ePPMs] = selfCal(Ps, Vs)

% [Ps_m, R_output1, T_output1] = self_calib_1(Ps, A, repmat(centers, nchoosek(10,2), 1));

% for i = 1:10
%   visualize(Ps_m(3*i-2:3*i,:))
% end


%% Visualize

function visualize(P)
x = [0  1  1  0  0  0  1  1  0  0  1  1  1  1  0  0];
y = [0  0  1  1  0  0  0  1  1  0  0  0  1  1  1  1];
z = [0  0  0  0  0  1  1  1  1  1  1  0  0  1  1  0];

figure
[m,n] = size(x);
x4d = [x(:),y(:),z(:),ones(m*n,1)]';
x2d = P(1:3,:)*x4d;
x2 = zeros(m,n); y2 = zeros(m,n);
x2(:) = x2d(1,:);
y2(:) = x2d(2,:);
plot(x2,y2)
end