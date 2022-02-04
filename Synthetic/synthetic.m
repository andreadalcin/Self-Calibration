clc, clear, close all
format long;

% Load parameters for sampling f
load('Synthetic/Data/sturm.mat')

% Load camera parameters
% K - intrinsic parameters (gt)
% width - image width
% height - image height
load('Synthetic/Data/camera_canon.mat')

% Noise in synthetic images
sigma = 0;

% Number of funamental matrices to sample - must be even
num_F = 10;

% Number of bins in initial guess histogram
num_bins = 100;

% Initial guess - principal points
u0 = 3072 / 2;
v0 = 2048 / 2;


%% Step 1 - Sample fundamental matrices

% Fs = zeros(3,3,num_F);
% for i = 1:num_F
%     %     [F1, F2] = sampleFundamental(K, sigma);
%     %     Fs(:,:,2*i-1) = F1;
%     %     Fs(:,:,2*i) = F2;
%     yaw = rand(1) * 45;
%     roll = rand(1) * 45;
%     pitch = rand(1) * 45;
%     R = rotYPR(yaw, roll, pitch);
% 
%     P1 = [K, [10, 5, 20]'];
%     P2 = K * [R, [1, 1, 1]'];
%     F = fund(P1, P2);
%     Fs(:,:,i) = F;
% end
% Fs_array = Fs;

load('Synthetic/Data/F.mat')
Fs_array = zeros(3,3,nchoosek(5,2));
index = 1;
for i = 1:size(Fs,3)-1
    for j = i+1:size(Fs,3)
        Fs_array(:,:,index) = Fs(:,:,i,j);
        index = index + 1;
    end
end


%% Step 2 - Get an initial estimate of K - K0

focal_lengths = [];

for i = 1:size(Fs_array,3)
    F = Fs_array(:,:,i);
    focal_lengths = [focal_lengths; getInitialEstimateOfFocalLength(F, width, height, opening_angles, false)];
end

h = histogram(focal_lengths, num_bins);
[~, mu] = robustcov(focal_lengths);

f0 = mu;
K0 = real([f0, 0, u0; 0, f0, v0; 0, 0, 1]);
disp('Intrinics by Sturm')
disp(K0)


%% Step 3 - Self-calibrate

X0 = [K0(1,:) K0(2,2:3)];

Options = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20,'MaxFunctionEvaluations',1e6);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fs_array, X, '2'), X0, [], [], Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

disp('Intrinsics computed by Mendonca&Cipolla`s Self Calibration')
disp(K_SK)