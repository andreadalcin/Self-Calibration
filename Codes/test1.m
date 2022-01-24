clear, clc;

format long;

% Load data
load('Fs.mat')
load('K_Canon.mat');
K_approx = K;
disp(K_approx)

% Intrinsics in Vector Form
X0 = [K_approx(1,:) K_approx(2,2:3)];
X0 = X0 + 10^2 * randn(1,5);
disp(X0);

%% Self Calibration using Simplified Kruppa's Method to Find Intrinsics

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionKruppaSimplified3(Fs, X),X0,[],[],Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];
disp('Intrinsics computed by Simplified Kruppa`s Self Calibration')
disp(K_SK)