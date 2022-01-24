function [F1,F2] = sampleFundamental(K, sigma)
num_points = 5*1e2;
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
t = [0 0 -10];

uv_kh = zeros(2,num_points,num_objects,num_images);

for i = 1:num_images
    % yaw = min_yaw + (max_yaw - min_yaw) * rand(1);
    % pitch = min_pitch + (max_pitch - min_pitch) * rand(1);
    uv_kh(:,:,:,i) = proj2dPoints(points_k, K, R, t) + sigma * randn(2,num_points,num_objects);
end

% Compute F1 and F2 from the two images
m1 = uv_kh(:,:,1,1);
m2 = uv_kh(:,:,1,2);
F0 = fund_lin(m2,m1,[]);
F1 = fund_nonlin(F0,m2,m1);

% [F1, inliers_1] = estimateFundamentalMatrix(m1', m2', 'Method', 'MSAC', 'NumTrials',2000, 'DistanceThreshold',1e-4);

% disp("F_1 check: x'F_1x = 0")
% res = [uv_kh(:,:,1,2); ones(1,num_points)]' * F1 * [uv_kh(:,:,1,1); ones(1,num_points)];
% disp(sum(diag(res)))

m1 = uv_kh(:,:,2,1);
m2 = uv_kh(:,:,2,2);
F0 = fund_lin(m2,m1,[]);
F2 = fund_nonlin(F0,m2,m1);

% [F2, inliers_2] = estimateFundamentalMatrix(m1', m2', 'Method', 'MSAC', 'NumTrials',2000, 'DistanceThreshold',1e-4);

% disp("F_2 check: x'F_2x = 0")
% res = [uv_kh(:,:,2,2); ones(1,num_points)]' * F2 * [uv_kh(:,:,2,1); ones(1,num_points)];
% disp(sum(diag(res)))



% Essential matrix
% disp("E_1 check: eigs")
% E_1 = K' * F1 * K;
% disp(eig(E_1))
% disp("E_2 check: eigs")
% E_2 = K' * F2 * K;
% disp(eig(E_2))

% K_approx = K + 10^0 * randn(3,3);
% norm((K - K_approx), 'fro')
% E_1 = K_approx' * F_1 * K_approx;
% eigs(E_1)
% E_2 = K_approx' * F_2 * K_approx;
% eigs(E_2)
end
