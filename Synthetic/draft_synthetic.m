clc, clear, close all
format long;

addpath('Synthetic/prior/')

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

outlier_ratio = 0;


%% Step 1 - Sample fundamental matrices
num_cams = 10;
sigma = 0;
outlier_indices = [];

index = 1;
for i = 1:num_cams-1
    for j = i+1:num_cams
        if (rand > outlier_ratio)
            K_curr = K;
        else
            a = 1e1;
            b = 1e5;
            r = (b-a).*rand(1) + a;
            K_curr = [r 0 u0; 0 r v0; 0 0 1];
            outlier_indices = [outlier_indices; index; index + 1];
        end

        Fs(:,:,index) = sampleFundamental(K,sigma);
        Fs(:,:,index+1) = sampleFundamental(K,sigma);

        index = index + 2;
    end
end


%% Step 2 - Get an initial estimate of K - K0

focal_lengths = [];

for i = 1:size(Fs,3)
    F = Fs(:,:,i);
    f = getInitialEstimateOfFocalLength(F, width, height);
    focal_lengths = [focal_lengths; f];
end

h = histogram(focal_lengths, num_bins);
[maxcount, whichbin] = max(h.Values);



[~,mu,~,~] = robustcov(focal_lengths);

% f0 = (h.BinEdges(whichbin) + h.BinEdges(whichbin+1)) / 2;
f0 = mu;
K0 = real([f0, 0, u0; 0, f0, v0; 0, 0, 1]);
disp('Intrinics by Sturm')
disp(K0)


%% Step 2.1 - Filter F
diffs = [];

for i = 1:size(Fs,3)
    F = Fs(:,:,i);
    estimates = getInitialEstimateFromF0(F, f0, u0, v0);
    diffs = [diffs; abs(K0(1,1) - median(estimates))];
end

figure;
plot(diffs, 'r.');

%%
threshold = prctile(diffs,10);
outliers = diffs > threshold;

disp("All outliers removed?")
sum(outliers(outlier_indices)) == size(outlier_indices,1)
disp("Percent:")
sum(outliers(outlier_indices)) / size(outlier_indices,1)

Fs_no_outliers = Fs;
Fs_no_outliers(:,:,outliers) = [];
disp("Remaining Fs:")
size(Fs_no_outliers,3)

% diffs = [];
% for i = 1:size(F_arr_no_outliers,3)
%     F = F_arr_no_outliers(:,:,i);
%     estimates = getInitialEstimateOfFocalLength(F, width, height, opening_angles, false);
%     if ~isempty(estimates)
%         diffs = [diffs; abs(K0(1,1) - median(estimates))];
%     end
% end
%
% figure;
% plot(diffs, 'b.');


%% Step 3 - Self-calibrate

X0 = [K0(1,:) K0(2,2:3)];

Options = optimoptions('lsqnonlin','Display','off',...
    'Algorithm','levenberg-marquardt',...
    'TolFun', 1e-20,'TolX',1e-20,...
    'MaxFunctionEvaluations',1e6);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionMendoncaCipollaVanilla(Fs_no_outliers, X, '2'), X0, [], [], Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

disp('Intrinsics computed by Mendonca & Cipolla`s Self Calibration')
disp(K_SK)


%% Step 4 - Check compatibility of K_SK with Fs

% for i = 1:size(Fs,3)
%     E = K_SK' * Fs(:,:,i) * K_SK;
%     disp(eigs(E))
% end





