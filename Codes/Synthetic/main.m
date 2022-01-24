clc, clear, close all

% rng('default')
load('../K_Canon.mat')
% K = eye(3);

num_points = 1e2;
num_objects = 2;    % k
num_images = 2;     % h

yaw = 45;
pitch = 90;

% Generate 3D model
points_k = zeros(3,num_points,num_objects);
for i = 1:num_objects
    [X,Y,Z] = randSph(num_points);
    points_k(:,:,i) = [X; Y; Z];
end

% Generate images
R = rotYPR(yaw,pitch,0);
t = [-10 -10 20];

uv_kh = zeros(2,num_points,num_objects,num_images);

for i = 1:num_images
    % yaw = min_yaw + (max_yaw - min_yaw) * rand(1);
    % pitch = min_pitch + (max_pitch - min_pitch) * rand(1);

    uv_kh(:,:,:,i) = proj2dPoints(points_k, K, R, t);
end

% Compute F1 and F2 from the two images
[F_1, inliers_1] = estimateFundamentalMatrix(uv_kh(:,:,1,1)', uv_kh(:,:,1,2)', 'Method', 'MSAC', 'NumTrials',2000, 'DistanceThreshold',1e-4);
[F_2, inliers_2] = estimateFundamentalMatrix(uv_kh(:,:,2,1)', uv_kh(:,:,2,2)', 'Method', 'MSAC', 'NumTrials',2000, 'DistanceThreshold',1e-4);

res = [uv_kh(:,:,1,2); ones(1,100)]' * F_1 * [uv_kh(:,:,1,1); ones(1,100)];
disp("F_1 check: x'F_1x = 0")
disp(sum(diag(res)))

res = [uv_kh(:,:,2,2); ones(1,100)]' * F_2 * [uv_kh(:,:,2,1); ones(1,100)];
disp("F_2 check: x'F_2x = 0")
disp(sum(diag(res)))

% Essential matrix
disp("E_1 check: eigs")
E_1 = K' * F_1 * K;
disp(eigs(E_1,2))
disp("E_2 check: eigs")
E_2 = K' * F_2 * K;
disp(eigs(E_2,2))

% K_approx = K + 10^0 * randn(3,3);
% norm((K - K_approx), 'fro')
% E_1 = K_approx' * F_1 * K_approx;
% eigs(E_1)
% E_2 = K_approx' * F_2 * K_approx;
% eigs(E_2)
