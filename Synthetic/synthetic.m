clc, clear, close all
format long;

% Load parameters for sampling f
load('Synthetic/Data/sturm.mat')

% Load camera parameters
% K - intrinsic parameters (gt)
% width - image width
% height - image height
load('Synthetic/Data/camera_canon.mat')

% Number of bins in initial guess histogram
num_bins = 100;

% Initial guess - principal points
u0 = width / 2;
v0 = height / 2;


%% Step 1 - Sample fundamental matrices
num_cams = 10;
sigma = 0;

index = 1;
for i = 1:num_cams-1
    for j = i+1:num_cams
        R = rotYPR(rand(1) * 45, rand(1) * 45, rand(1) * 45);
        P1 = [K, [10, 5, 20]'];
        P2 = [K * R, [1, 1, 1]'];
        Fs(:,:,index) = fund(P1, P2);

        R = rotYPR(rand(1) * 45, rand(1) * 45, rand(1) * 45);
        P1 = [K, [10, 5, 20]'];
        P2 = [K * R, [1, 1, 1]'];
        Fs(:,:,index+1) = fund(P1, P2);
        %         [F1, F2] = sampleFundamental(K,sigma);
        %         Fs(:,:,index) = F1;
        %         Fs(:,:,index+1) = F2;

        index = index + 2;
    end
end


%% Step 2 - Get an initial estimate of K - K0

focal_lengths = [];

for i = 1:size(Fs,3)
    F = Fs(:,:,i);
    focal_lengths = [focal_lengths; getInitialEstimateOfFocalLength(F, width, height, opening_angles, false)];
end

h = histogram(focal_lengths, num_bins);
[~,mu,~,outliers] = robustcov(focal_lengths);

f0 = mu;
K0 = real([f0, 0, u0; 0, f0, v0; 0, 0, 1]);
disp('Intrinics by Sturm')
disp(K0)


%% Step 3 - Self-calibrate

X0 = [K0(1,:) K0(2,2:3)];

Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20,'MaxFunctionEvaluations',1e6);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fs, X, '1'), X0, [], [], Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

disp('Intrinsics computed by Mendonca & Cipolla`s Self Calibration')
disp(K_SK)


%% Step 4 - Check compatibility of K_SK with Fs
residual = 0;
for i = 1:size(Fs,3)
    E = K_SK' * Fs(:,:,i) * K_SK;
    M = eigs(E,2);
    residual = residual + 1 - real(M(2)) / real(M(1));
end

disp("Average residual")
disp(residual / size(Fs,3))




