clear, clc;

format long;

%% Load data


% load('./Data/F12.mat');
% Fs(:,:,1) = F1;
% Fs(:,:,2) = F2;
% 
% load('./Data/F14.mat');
% Fs(:,:,3) = F1;
% Fs(:,:,4) = F2;
% 
% load('./Data/F15.mat');
% Fs(:,:,5) = F1;
% Fs(:,:,6) = F2;
% 
% load('./Data/F23.mat');
% Fs(:,:,7) = F1;
% Fs(:,:,8) = F2;
% 
% load('./Data/F25.mat');
% Fs(:,:,9) = F1;
% Fs(:,:,10) = F2;

Fs = double(zeros(3,3,5));

load('./Data1/F13.mat')
Fs(:,:,1) = F1;
Fs(:,:,2) = F2;
load('./Data1/F23.mat')
Fs(:,:,3) = F1;
Fs(:,:,4) = F2;
load('./Data1/F12.mat')
Fs(:,:,5) = F1;
% Fs(:,:,6) = F2;

load('K_Canon.mat');
K_approx = K;
disp(K_approx)

% Intrinsics in Vector Form
X0 = [K_approx(1,:) K_approx(2,2:3)];
X0 = X0 + 10^2 * randn(1,5);
disp(X0);

%% Self Calibration using Simplified Kruppa's Method to Find Intrinsics

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt','TolX', 1e-10);

% Compute Intrinsics by Optimization(Method 1)
K_MC1 = lsqnonlin(@(X) costFunctionMendoncaCipolla2(Fs, X, '1'),X0,[],[],Options);

% Intrinsics in Matrix Form
K_MC1 = [K_MC1(1) K_MC1(2) K_MC1(3); 0 K_MC1(4) K_MC1(5); 0 0 1];
disp('Intrinsics computed by Mendonca & Cipolla`s Self Calibration(Method 1)')
disp(K_MC1)

% Compute Intrinsics by Optimization(Method 2)
K_MC2 = lsqnonlin(@(X) costFunctionMendoncaCipolla2(Fs, X, '2'),X0,[],[],Options);

% Intrinsics in Matrix Form
K_MC2 = [K_MC2(1) K_MC2(2) K_MC2(3); 0 K_MC2(4) K_MC2(5); 0 0 1];
disp('Intrinsics computed by Mendonca & Cipolla`s Self Calibration(Method 2)')
disp(K_MC2)