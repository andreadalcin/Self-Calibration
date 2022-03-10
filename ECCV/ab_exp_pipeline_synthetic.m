clc, clear, close all
format longG;

addpath('Synthetic/')
addpath('Synthetic/prior/')

% Load camera parameters
% K - intrinsic parameters (gt)
% width - image width
% height - image height
load('Synthetic/Data/camera_canon.mat')

% Parameters
num_cams = 4;
sigma = 0;
num_samples = 1;
motions = 1;
min_outlier_f = 1e1;
max_outlier_f = 1e5;
outlier_ratio = 0.0;

Fx = [];
Fy = [];
U = [];
V = [];

for i = 1:num_samples
    i
    [fx, fy, u, v] = run_test(K, width, height, motions, num_cams, sigma, ...
        min_outlier_f, max_outlier_f, outlier_ratio);
    Fx = [Fx; fx];
    Fy = [Fy; fy];
    U = [U; u];
    V = [V; v];
end

median(Fx)
median(Fy)
median(U)
median(V)


function [fx, fy, u, v] = run_test(K, width, height, motions, ...
    num_cams, sigma, min_outlier_f, max_outlier_f, outlier_ratio)

load('Dataset/synthetic/test-7-3.mat')

% Sample fundamental matrices
% Fs = [];
% for i = 1:num_cams-1
%     for j = i+1:num_cams
%         for k = 1:motions
%             if (rand > outlier_ratio)
%                 K_curr = K;
%             else
%                 r = (max_outlier_f - min_outlier_f) .* rand(1) + min_outlier_f;
%                 K_curr = [r 0 width / 2; 0 r height / 2; 0 0 1];
%             end
% 
%             Fs(:,:,size(Fs,3) + 1) = sampleFundamental(K_curr, sigma);
%         end
%     end
% end
% Fs(:,:,1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the initial focal length
[mu_f0, sigma_f0] = sturm2(Fs, width, height);
fprintf("mu_f0: %f, sigma_f0: %f\n", mu_f0, sigma_f0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = ones(size(Fs,3),1);
weights = weights ./ sum(weights);

Fx = [];
Fy = [];
views = [];

% residual = inf;
% res_fx = inf;
% res_fy = inf;

num_trials = 1000;
for i = 1:num_trials
    [fx, fy, v] = optimizeFocalLength(Fs, weights, mu_f0, sigma_f0 * 0.1, width, height);
    Fx(i) = fx;
    Fy(i) = fy;
    views(i,:) = v;

%     K_curr = [fx 0 width / 2; 0 fy height / 2; 0 0 1];
%     curr_residual = costFunctionMendoncaCipollaVanilla(Fs, K_curr, '2');
%     if residual > curr_residual
%         residual = curr_residual;
%         res_fx = fx;
%         res_fy = fy;
%     end
end

fx = kernel_voting(Fx', 0.01);
fy = kernel_voting(Fy', 0.01);

% Find k-NN
kNN = num_trials * 0.01;
I = knnsearch([Fx' Fy'], [fx fy], 'K', kNN);
v = views(I,:);
[ii,~,kk] = unique(v(:));
v = ii(histc(kk,1:numel(ii)) > 1);

% Refinement with KNN
[fx, fy, u, v] = optimizePrincipalPoints(Fs(:,:,v), fx, fy, width, height);

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
    'MaxIterations',1e2,...
    'MaxFunctionEvaluations',1e6,...
    'TolFun', 1e-20,...
    'TolX',1e-20);

% Randomize focal length
X0 = [K0(1,:) K0(2,2:3)];

loss = @(X) costFunctionMendoncaCipollaVanilla(Fs, X, '2');

K_SK = lsqnonlin(loss, X0, [], [], Options);
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

fx = K_SK(1,1);
fy = K_SK(2,2);
u = K_SK(1,3);
v = K_SK(2,3);

end
