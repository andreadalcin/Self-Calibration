function M = sampleMatches(K, sigma, num_bodies)

num_points = 200;
num_objects = 3;    % k
num_images = 2;     % h

% Generate 3D model
points_k = zeros(3,num_points,num_objects);
for i = 1:num_objects
    [X,Y,Z] = randSph(num_points);
    points_k(:,:,i) = [X; Y; Z];
end

% Generate images
yaw = 90;
pitch = 90;
roll = 0;
R = rotYPR(yaw,pitch,roll);
t = [0 0 -10];

uv_kh = zeros(2,num_points,num_objects,num_images);

for i = 1:num_images
    uv_kh(:,:,:,i) = proj2dPoints(points_k, K, R, t) + sigma * randn(2,num_points,num_objects);
end

% Compute F1 and F2 from the two images
src_1 = uv_kh(:,:,1,1);
dst_1 = uv_kh(:,:,1,2);
src_2 = uv_kh(:,:,2,1);
dst_2 = uv_kh(:,:,2,2);
src_3 = uv_kh(:,:,3,1);
dst_3 = uv_kh(:,:,3,2);

M = zeros(3 * num_points,4);
M(1:num_points,1:2) = src_1';
M(1:num_points,3:4) = dst_1';
M(num_points+1:2*num_points,1:2) = src_2';
M(num_points+1:2*num_points,3:4) = dst_2';
M(2*num_points+1:3*num_points,1:2) = src_3';
M(2*num_points+1:3*num_points,3:4) = dst_3';

M = M(1:num_points*num_bodies,:,:);


% F0 = fund_lin(m2,m1,[]);
% F1 = fund_nonlin(F0,m2,m1);
% [F1, ~] = estimateFundamentalMatrix(m1', m2');

% disp("F_1 check: x'F_1x = 0")
% res = [uv_kh(:,:,1,2); ones(1,num_points)]' * F1 * [uv_kh(:,:,1,1); ones(1,num_points)];
% disp(sum(diag(res)))

% m1 = uv_kh(:,:,2,1);
% m2 = uv_kh(:,:,2,2);
% F0 = fund_lin(m2,m1,[]);
% F2 = fund_nonlin(F0,m2,m1);
% [F2, ~] = estimateFundamentalMatrix(m1', m2', 'Method', 'MSAC', 'NumTrials', 2000, 'DistanceThreshold', 1e-4);

% disp("F_2 check: x'F_2x = 0")
% res = [uv_kh(:,:,2,2); ones(1,num_points)]' * F2 * [uv_kh(:,:,2,1); ones(1,num_points)];
% disp(sum(diag(res)))

end
