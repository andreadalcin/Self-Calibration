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

num_cams = 10;
sigma = 0;

min_outlier_f = 1e1;
max_outlier_f = 1e5;
outlier_ratio = 0.5;


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

%% Step 1.1 - Load a real dataset

load("ECCV/Data1/F_fountain.mat", 'Fs')
num_cameras = size(Fs,3);
width = 3072;
height = 2048;
Fs_f = [];
index = 1;
for i = 1:num_cameras-1
    for j = i+1:num_cameras
        if Fs(:,:,i,j,1) ~= zeros(3,3)
            Fs_f(:,:,index) = Fs(:,:,i,j,1);
            index = index + 1;
            % Fs_f(:,:,index) = Fs(:,:,i,j,2);
            % index = index + 1;
        end
    end
end
Fs_sb = Fs_f;
Fs_mb = Fs_f;
outliers_sb = [];
outliers_mb = [];


%% Step 2 - Evaluate

pct0 = 20;

disp("Single body")
evaluate(Fs_sb, outliers_sb, K, width, height, num_bins, pct0);



function [d_f0, d_f, d_uv] = evaluate(Fs, outliers, K, width, height, num_bins, pct0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [];
for i = 1:size(Fs,3)
    f0 = getInitialEstimateOfFocalLength(Fs(:,:,i), width, height);
    f = [f; f0];
end

figure;
h = histogram(f, num_bins);

% Filter focal lengths not in the peak bin
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);
f(or(f < th_low, f > th_high)) = [];

f0 = median(f);
fprintf("f0: %d\n", f0)
K0 = [f0, 0, width / 2; 0, f0, height / 2; 0, 0, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residual = [];
for i = 1:size(Fs,3)
    f = getInitialEstimateFromF0(Fs(:,:,i), f0, width, height);
    if isempty(f)
        residual = [residual; 1e5];
    else
        residual = [residual; abs(f0 - f)];
    end
end

Fo = [];
pct = pct0;
while (size(Fo,3) <= 6)
    th_low = prctile(residual,pct);
    Fo = Fs(:,:,residual <= th_low);
    pct = pct + 5;
end

figure;
hold on
plot(residual, 'r.');
yline(th_low);
hold off

disp("Outlier removal percentage:")
o = residual > th_low;
sum(o(outliers)) / size(outliers,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options = optimoptions('lsqnonlin','Display','','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20,'MaxFunctionEvaluations',1e6);

K0 = real([f0, 0, width / 2; 0, f0, height / 2; 0, 0, 1]);
X0 = [K0(1,:) K0(2,2:3)];
K_SK = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fo, X, '1'), X0, [], [], Options);
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

disp('Intrinsics - Mendonca&Cipolla')
disp(K_SK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_f0 = (abs(K(1,1) - f0) + abs(K(2,2) - f0)) / 2;
d_f = (abs(K(1,1) - K_SK(1,1)) + abs(K(2,2) - K_SK(2,2))) / 2;
d_uv = norm([K(1,3) - K_SK(1,3), K(2,3) - K_SK(2,3)]);
disp("Number of Fundamental matrices")
disp(size(Fs,3))
disp("Initial focal length - error")
disp(d_f0)
disp("Focal length - error")
disp(d_f)
disp("Principal point - error")
disp(d_uv)

end
