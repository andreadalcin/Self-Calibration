clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

load('Dataset/Fs.mat')

width = 4272;
height = 2848;

A = [5e3 0 width/2; 0 5e3 height / 2; 0 0 1];


%% Self Calibration using Dual Absolute Quadric Method to Find Intrinsics

P = vgg_P_from_F(Fs(:,:,1,3,2));

X0 = [A(1,:) A(2,:) A(3,:)];

% Compute Normal to Plane at Infinity
Normal = NormaltoPlaneAtInfinity(P, A);

% Compute Homography of Plane at Infinity
H = homographyPlaneAtInfinity(Fs(:,:,:,:,1), Normal);

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter','TolFun', 1e-16,'TolX',1e-16);

% Compute Intrinsics by Optimization
K_DAQ = lsqnonlin(@(X) costFunctionDualAbsoluteQuadric(H, X), X0, [], [], Options);

% Intrinsics in Matrix Form
K_DAQ = [K_DAQ(1) K_DAQ(2) K_DAQ(3); K_DAQ(4) K_DAQ(5) K_DAQ(6); K_DAQ(7) K_DAQ(8) K_DAQ(9)];

% Regain the Scale by setting the Last element to 1
K_DAQ = (1/K_DAQ(3,3))*K_DAQ;
disp('Intrinsics computed by DAQ Self Calibration')
disp(K_DAQ)