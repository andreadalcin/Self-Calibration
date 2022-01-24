clc, clear, close all
format long;
addpath('../')
run("ComputerVisionToolkit/cvt_setup.m");

% rng('default')

% Load data
load('../data.mat')
K = A;

sigma = 1e-2;
num_iter = 10;

Fs = zeros(3,3,2*num_iter);

for i = 1:num_iter
    [F1, F2] = sampleFundamental(K, sigma);
    Fs(:,:,2*i-1) = F1;
    Fs(:,:,2*i) = F2;
end

% Intrinsics in Vector Form
K_approx = K;
K_approx = K_approx + rand(3,3);
K_approx(2,1) = 0;
K_approx(1,2) = 0;
K_approx(3,1) = 0;
K_approx(3,2) = 0;
K_approx = K_approx ./ K_approx(3,3);
X0 = [K_approx(1,:) K_approx(2,2:3)];

disp('K:')
disp(K)
disp('K_approx (initial guess):')
disp(K_approx)

% Self Calibration using Simplified Kruppa's Method to Find Intrinsics

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20,'MaxFunctionEvaluations',1e6);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionKruppaSimplified(Fs, X),X0,[],[],Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

K_SK(1,2) = 0;

disp('Intrinsics computed by Simplified Kruppa`s Self Calibration')
disp(K_SK)
disp('Frobenius:')
disp(norm(K_SK - K, "fro"))
