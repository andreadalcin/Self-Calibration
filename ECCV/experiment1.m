clc, clear, close all
format longG;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

dataset = "fountain";
% dataset = "herzjesu";
% dataset = "castle";
% dataset = "01_amiibo_static";
% dataset = "02_amiibo_motion";
% dataset = "04_amiibo_motion";
% dataset = "06_amiibo_static";
% dataset = "10_amiibo_motion";
% dataset = "12_amiibo_motion";
% dataset = "KITTI_Seq005";
% dataset = "KITTI_Seq113";

switch dataset
    case "fountain"
        load('Dataset/Sturm/fountain.mat')
        width = 3072;
        height = 2048;
        fx_gt = 2759.48;
        fy_gt = 2764.16;
        u_gt = 1520.69;
        v_gt = 1006.81;
        Fs = pre_process(Fs);
    case "01_amiibo_static"
        load('Dataset/Sturm/01_amiibo_static.mat')
        width = cameraParams.ImageSize(2);
        height = cameraParams.ImageSize(1);
        fx_gt = cameraParams.Intrinsics.FocalLength(1);
        fy_gt = cameraParams.Intrinsics.FocalLength(2);
        u_gt = cameraParams.Intrinsics.PrincipalPoint(1);
        v_gt = cameraParams.Intrinsics.PrincipalPoint(2);
        Fs = pre_process(Fs);
end


load('Dataset/optimization/10_amiibo_motion.mat')
load('Dataset/params.mat')

% Load parameters

matches = [];

kNN_support = 10;
num_cameras = size(Fs,3);

% Experiments
sample_size = 1;
num_trials = 5000;

Fx = []; % Fx = [num_cameras, SB, MB]
Fy = []; % Fy = [num_cameras, SB, MB]

for num = num_cameras:num_cameras
    Fx_sb = [];
    Fy_sb = [];
    U_sb = [];
    V_sb = [];

    Fx_mb = [];
    Fy_mb = [];
    U_mb = [];
    V_mb = [];

    for i = 1:sample_size
        cameras = sort(datasample(1:size(Fs,3), num, 'Replace', false))

        [fx, fy, u, v] = selfCalibrate(false, cameras, Fs, ...,
            width, height, num_trials, kNN_support, [], []);
        Fx_sb = [Fx_sb; fx];
        Fy_sb = [Fy_sb; fy];
        U_sb = [U_sb; u];
        V_sb = [V_sb; v];

        [fx, fy, u, v] = selfCalibrate(true, cameras, Fs, ...
            width, height, num_trials, kNN_support, [], []);
        Fx_mb = [Fx_mb; fx];
        Fy_mb = [Fy_mb; fy];
        U_mb = [U_mb; u];
        V_mb = [V_mb; v];
    end

    % Single body
    fx = median(Fx_sb);
    fy = median(Fy_sb);
    u = median(U_sb);
    v = median(V_sb);
    err_f_sb = 0.5 * 100 * (abs((fx - fx_gt)/fx_gt) + abs((fy - fy_gt)/fy_gt));
    err_pp_sb = 0.5 * 100 * (abs((u - u_gt)/u_gt) + abs((v - v_gt)/v_gt));

    % Multi body
    fx = median(Fx_mb);
    fy = median(Fy_mb);
    u = median(U_mb);
    v = median(V_mb);
    err_f_mb = 0.5 * 100 * (abs((fx - fx_gt)/fx_gt) + abs((fy - fy_gt)/fy_gt));
    err_pp_mb = 0.5 * 100 * (abs((u - u_gt)/u_gt) + abs((v - v_gt)/v_gt));

    disp("df - SB")
    disp(err_f_sb)
    disp("df - MB")
    disp(err_f_mb)
    disp("dpp - SB")
    disp(err_pp_sb)
    disp("dpp - MB")
    disp(err_pp_mb)

    Fx = [Fx; [num_cameras, ...
        abs(cameraParams.Intrinsics.FocalLength(1) - median(Fx_sb)), ...
        var(cameraParams.Intrinsics.FocalLength(1) - Fx_sb), ...
        abs(cameraParams.Intrinsics.FocalLength(1) - median(Fx_mb)), ...
        var(cameraParams.Intrinsics.FocalLength(1) - Fx_mb)]];

    Fy = [Fy; [num_cameras, ...
        abs(cameraParams.Intrinsics.FocalLength(2) - median(Fy_sb)), ...
        var(cameraParams.Intrinsics.FocalLength(2) - Fy_sb), ...
        abs(cameraParams.Intrinsics.FocalLength(2) - median(Fy_mb)), ...
        var(cameraParams.Intrinsics.FocalLength(2) - Fy_mb)]];
end


%%
figure
title("FX")
hold on
plot(Fx(:,1), Fx(:,2), 'r');
plot(Fx(:,1), Fx(:,4), 'b');
hold off

figure
title("FY")
hold on
plot(Fy(:,1), Fy(:,2), 'r');
plot(Fy(:,1), Fy(:,5), 'b');
hold off


function [fx, fy, u, v] = selfCalibrate(useMultibody, cameras, Fs, ...
    width, height, numTrials, knnSupport, matches_1, matches_2)

num_cameras = size(cameras,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fundamental matrices and priors
Fo = [];
priors = [];
for a = 1:num_cameras-1
    for b = a+1:num_cameras
        i = cameras(a);
        j = cameras(b);

        if Fs(:,:,i,j,1) ~= zeros(3,3)
            Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,1) / norm(Fs(:,:,i,j,1));
        end

        if useMultibody
            if Fs(:,:,i,j,2) ~= zeros(3,3)
                Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,2) / norm(Fs(:,:,i,j,2));
            end
        end
    end
end
Fo(:,:,1) = [];

% priors = priors ./ sum(priors);
priors = ones(size(Fo,3),1);
priors = priors ./ sum(priors);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the initial focal length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mu_f0, sigma_f0] = sturm2(Fo, width, height);
fprintf("mu_f0: %f, sigma_f0: %f\n", mu_f0, sigma_f0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = ones(size(Fo,3),1);
weights = weights ./ sum(weights);

Fx = [];
Fy = [];
views = [];

% residuals = zeros(numTrials,1);
% K_residual = zeros(3,3,numTrials);

parfor i = 1:numTrials
    [fx, fy, v] = optimizeFocalLength(Fo, weights, mu_f0, sigma_f0 * 0.1, width, height);
    Fx(i) = fx;
    Fy(i) = fy;
    views(i,:) = v;

    % K_residual(:,:,i) = [fx, 0, width/2; 0 fy height/2; 0 0 1];
    % residuals(i) = sum(costFunctionMendoncaCipollaResidual(Fo, K_residual(:,:,i), priors, '2'));
end

% [~,I] = min(residuals);
% K_residual = K_residual(:,:,I);
%
% fx = K_residual(1,1)
% fy = K_residual(2,2)

fx = kernel_voting(Fx', 0.05);
fy = kernel_voting(Fy', 0.05);

% Find k-NN
kNN = numTrials * 0.01;
I = knnsearch([Fx' Fy'], [fx fy], 'K', kNN);
v = views(I,:);
[ii,~,kk] = unique(v(:));
v = ii(histc(kk,1:numel(ii)) > knnSupport);

% Refinement with KNN
[fx, fy, u, v] = optimizePrincipalPoints(Fo(:,:,v), fx, fy, width, height);

end


function f = kernel_voting(f, bandwidth_factor)
f(f < 100 | f > 1e5) = [];
bandwidth = median(f) * bandwidth_factor;
pdSix = fitdist(f,'Kernel','Width',bandwidth);
x = min(f):.1:max(f);
ySix = pdf(pdSix,x);
[~,I] = max(ySix);
f = x(I);
figure
plot(x,ySix)
end


function [fx,fy,views] = optimizeFocalLength(Fs, weights, mu_f0, sigma_f0, width, height)

% Randomize focal length and principal point
f0 = 0;
u0 = 0;
v0 = 0;

while (f0 <= 0)
    f0 = normrnd(mu_f0, sigma_f0);
end
while (u0 <= 0)
    u0 = normrnd(width / 2, width / 6);
end
while (v0 <= 0)
    v0 = normrnd(height / 2, height / 6);
end

K0 = [f0, 0, u0; 0, f0, v0; 0, 0, 1];

% Optimization
Options = optimoptions('lsqnonlin','Display','off', ...
    'Algorithm','levenberg-marquardt', ...
    'StepTolerance',1e-20,...
    'FunctionTolerance',1e-20,...
    'MaxIterations',1e2,...
    'MaxFunctionEvaluations',1e6,...
    'TolFun', 1e-20,...
    'TolX',1e-20);

% Randomize focal length
X0 = [K0(1,1) K0(2,2)];

views = datasample(1:size(Fs,3),3,'Replace',false,'Weights',weights);
subset = Fs(:,:,views);

loss = @(X) costFunctionMendoncaCipollaFocalOnly(subset, X, u0, v0, '2');

K_SK = lsqnonlin(loss, X0, [], [], Options);
K_SK = [K_SK(1) 0 width/2; 0 K_SK(2) height/2; 0 0 1];

fx = K_SK(1,1);
fy = K_SK(2,2);

end


function [fx,fy,u,v] = optimizePrincipalPoints(Fs, fx, fy, width, height)

u0 = width / 2;
v0 = height / 2;
K0 = [fx, 0, u0; 0, fy, v0; 0, 0, 1];

% Optimization
Options = optimoptions('lsqnonlin','Display','off', ...
    'Algorithm','levenberg-marquardt', ...
    'StepTolerance',1e-20,...
    'FunctionTolerance',1e-20,...
    'MaxIterations',10,...
    'MaxFunctionEvaluations',1e6,...
    'TolFun', 1e-20,...
    'TolX',1e-20);

% Randomize focal length
X0 = [K0(1,1) K0(2,2) K0(1,3), K0(2,3)];

loss = @(X) costFunctionMendoncaCipolla(Fs, X, '2');

K_SK = lsqnonlin(loss, X0, [], [], Options);
K_SK = [K_SK(1) 0 K_SK(3); 0 K_SK(2) K_SK(4); 0 0 1];

fx = K_SK(1,1);
fy = K_SK(2,2);
u = K_SK(1,3);
v = K_SK(2,3);

end

function [Fo, weights] = pre_process(Fs)
Fo = [];
weights = [];
num_cameras = size(Fs,3);
for i = 1:num_cameras-1
    for j = i+1:num_cameras
        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,1) / norm(Fs(:,:,i,j,1));
    end
end
Fo(:,:,1) = [];
end

function [Fo, weights] = pre_process_kitti(Fs)
Fo = [];
weights = [];
for i = 1:size(Fs,3)-1
    for j = i+1:size(Fs,3)
        if Fs(:,:,i,j) ~= zeros(3,3)
            Fo(:,:,size(Fo,3) + 1) = Fs(:,:,i,j);
        end
    end
end
Fo(:,:,1) = [];
end