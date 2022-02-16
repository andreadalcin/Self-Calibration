clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')

% Load camera parameters
% K - intrinsic parameters (gt)
% width - image width
% height - image height
load('Synthetic/Data/camera_canon.mat')

% Parameters
num_bins = 50;

translation_scale = 1e3;
rotation_scale = 45;

num_cams = 5;
sigma = 0;

min_outlier_f = 1e1;
max_outlier_f = 1e5;
outlier_ratio = 0.6;


%% Step 1 - Sample fundamental matrices
Fs_sb = zeros(3,3,nchoosek(num_cams,2));
Fs_mb = zeros(3,3,2*nchoosek(num_cams,2));

outliers_sb = [];
outliers_mb = [];

i_sb = 1;
i_mb = 1;
for i = 1:num_cams-1
    for j = i+1:num_cams
        % Sample first fundamental matrix
        if (rand > outlier_ratio)
            K_curr = K;
        else
            r = (max_outlier_f - min_outlier_f) .* rand(1) + min_outlier_f;
            K_curr = [r 0 width / 2; 0 r height / 2; 0 0 1];
            outliers_sb = [outliers_sb; i_sb];
            outliers_mb = [outliers_mb; i_mb];
        end

        P1 = K_curr * [eye(3,3), zeros(3,1)];
        R = rotYPR(rand(1) * rotation_scale, rand(1) * rotation_scale, rand(1) * rotation_scale);
        t = translation_scale * rand(3,1);
        P2 = K_curr * [R, t];
        F = fund(P1, P2);
        % F = sampleFundamental(K_curr, sigma);
        Fs_sb(:,:,i_sb) = F;
        Fs_mb(:,:,i_mb) = F;
        i_sb = i_sb + 1;
        i_mb = i_mb + 1;

        % Sample second fundamental matrix (MB only)
        if (rand > outlier_ratio)
            K_curr = K;
        else
            r = (max_outlier_f - min_outlier_f) .* rand(1) + min_outlier_f;
            K_curr = [r 0 width / 2; 0 r height / 2; 0 0 1];
            outliers_mb = [outliers_mb; i_mb];
        end

        P1 = K_curr * [eye(3,3), zeros(3,1)];
        R = rotYPR(rand(1) * rotation_scale, rand(1) * rotation_scale, rand(1) * rotation_scale);
        t = translation_scale * rand(3,1);
        P2 = K_curr * [R, t];
        F = fund(P1, P2);
        % F = sampleFundamental(K_curr, sigma);
        Fs_mb(:,:,i_mb) = F;
        i_mb = i_mb + 1;
    end
end


%% Step 2 - Estimate focal length from fundamental matrices

% Single body
f_sb = [];
for i = 1:size(Fs_sb,3)
    f0 = getInitialEstimateOfFocalLength(Fs_sb(:,:,i), width, height);
    f_sb = [f_sb; f0];
end

h = histogram(f_sb, num_bins);
[maxcount, whichbin] = max(h.Values);
f0_sb = (h.BinEdges(whichbin) + h.BinEdges(whichbin+1)) / 2;
% [~,f0_sb,~,~] = robustcov(f_sb);
K0_sb = [f0_sb, 0, width/2; 0, f0_sb, height/2; 0, 0, 1];
disp('f0 - Single body')
disp(f0_sb)

% Multi-body
f_mb = [];
for i = 1:size(Fs_mb,3)
    f0 = getInitialEstimateOfFocalLength(Fs_mb(:,:,i), width, height);
    f_mb = [f_mb; f0];
end

h = histogram(f_mb, num_bins);
[maxcount, whichbin] = max(h.Values);
f0_mb = (h.BinEdges(whichbin) + h.BinEdges(whichbin+1)) / 2;
% [~,f0_mb,~,~] = robustcov(f_mb);
K0_sb = [f0_mb, 0, width/2; 0, f0_mb, height/2; 0, 0, 1];
disp('f0 - Multi body')
disp(f0_mb)


%% Step 4 - Remove outliers
residual_sb = [];
for i = 1:size(Fs_sb,3)
    f0 = getInitialEstimateFromF0(Fs_sb(:,:,i), f0_sb, width, height);
    if isempty(f0)
        residual_sb = [residual_sb; 1e5];
    else
        residual_sb = [residual_sb; abs(f0_sb - f0)];
    end
end

Fo_sb = [];
pct = 5;
while (size(Fo_sb,3) <= 6)
    th_low = prctile(residual_sb,pct);
    Fo_sb = Fs_sb(:,:,residual_sb <= th_low);
    pct = pct + 5;
end

figure;
hold on
plot(residual_sb, 'r.');
yline(th_low);
hold off

disp("SB - Outlier removal percentage:")
outliers = residual_sb > th_low;
sum(outliers(outliers_sb)) / size(outliers_sb,1)


residual_mb = [];
for i = 1:size(Fs_mb,3)
    f0 = getInitialEstimateFromF0(Fs_mb(:,:,i), f0_mb, width, height);
    if isempty(f0)
        residual_mb = [residual_mb; 1e9];
    else
        residual_mb = [residual_mb; abs(f0_mb - f0)];
    end
end

Fo_mb = [];
pct = 5;
while (size(Fo_mb,3) <= 6)
    th_low = prctile(residual_mb,pct);
    Fo_mb = Fs_mb(:,:,residual_mb <= th_low);
    pct = pct + 5;
end

figure;
hold on
plot(residual_mb, 'b.');
yline(th_low);
hold off

disp("MB - Outlier removal percentage:")
outliers = residual_mb > th_low;
sum(outliers(outliers_mb)) / size(outliers_mb,1)


%% Step 3 - Self-calibration - Single body

Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20,'MaxFunctionEvaluations',1e6);

K0 = real([f0_sb, 0, width / 2; 0, f0_sb, height / 2; 0, 0, 1]);
X0 = [K0(1,:) K0(2,2:3)];
K_sb = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fo_sb, X, '1'), X0, [], [], Options);
K_sb = [K_sb(1) K_sb(2) K_sb(3); 0 K_sb(4) K_sb(5); 0 0 1];

disp('Intrinsics - Single body')
disp(K_sb)

K0 = real([f0_mb, 0, width / 2; 0, f0_mb, height / 2; 0, 0, 1]);
X0 = [K0(1,:) K0(2,2:3)];
K_mb = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fo_mb, X, '1'), X0, [], [], Options);
K_mb = [K_mb(1) K_mb(2) K_mb(3); 0 K_mb(4) K_mb(5); 0 0 1];
disp('Intrinsics - Multi body')
disp(K_mb)


%% Step 4 - Metrics

disp("Number of cameras")
disp(num_cams)

% Single body
d_f0 = (abs(K(1,1) - f0_sb) + abs(K(2,2) - f0_sb)) / 2;
d_f = (abs(K(1,1) - K_sb(1,1)) + abs(K(2,2) - K_sb(2,2))) / 2;
d_uv = norm([K(1,3) - K_sb(1,3), K(2,3) - K_sb(2,3)]);
disp("Single body - Number of Fundamental matrices")
disp(size(Fs_sb,3))
disp("Single body - Initial focal length - error")
disp(d_f0)
disp("Single body - Focal length - error")
disp(d_f)
disp("Single body - Principal point - error")
disp(d_uv)

d_f0 = (abs(K(1,1) - f0_mb) + abs(K(2,2) - f0_mb)) / 2;
d_f = (abs(K(1,1) - K_mb(1,1)) + abs(K(2,2) - K_mb(2,2))) / 2;
d_uv = norm([K(1,3) - K_mb(1,3), K(2,3) - K_mb(2,3)]);
disp("Multi body - Number of Fundamental matrices")
disp(size(Fs_mb,3))
disp("Multi body - Initial focal length - error")
disp(d_f0)
disp("Multi body - Focal length - error")
disp(d_f)
disp("Multi body - Principal point - error")
disp(d_uv)
