clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

load('Dataset/fountain.mat')
load('Dataset/params.mat')

% Load parameters
fx_gt = cameraParams.Intrinsics.FocalLength(1);
fy_gt = cameraParams.Intrinsics.FocalLength(2);
u_gt = cameraParams.Intrinsics.PrincipalPoint(1);
v_gt = cameraParams.Intrinsics.PrincipalPoint(2);

% width = cameraParams.ImageSize(2);
% height = cameraParams.ImageSize(1);
matches = [];
width = 3072;
height = 2048;

num_cameras = size(Fs,3);

% Experiments
sample_size = 1;
num_trials = 100;

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
            matches, ...
            width, height, num_trials);
        Fx_sb = [Fx_sb; fx];
        Fy_sb = [Fy_sb; fy];
        U_sb = [U_sb; u];
        V_sb = [V_sb; v];
        close all

        [fx, fy, u, v] = selfCalibrate(true, cameras, Fs, ...
            matches, ...
            width, height, num_trials);
        Fx_mb = [Fx_mb; fx];
        Fy_mb = [Fy_mb; fy];
        U_mb = [U_mb; u];
        V_mb = [V_mb; v];
        close all
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


function [fx, fy, u, v] = selfCalibrate(useMultibody, cameras, Fs, matches, width, height, numTrials)

num_cameras = size(cameras,2);
p = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fundamental matrices and priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For KITTI
% Fo = Fs;
% Fo = [];
% for i = 1:size(Fs,3)-1
%     for j = i+1:size(Fs,3)
%         if Fs(:,:,i,j) ~= zeros(3,3)
%             Fo(:,:,size(Fo,3) + 1) = Fs(:,:,i,j);
%         end
%     end
% end
% Fo(:,:,1) = [];

% For others
Fo = [];
weights = [];

for a = 1:num_cameras-1
    for b = a+1:num_cameras
        i = cameras(a);
        j = cameras(b);

        if rand > p
            o1 = 1;
            o2 = 2;
        else
            o1 = 2;
            o2 = 1;
        end

        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o1) / norm(Fs(:,:,i,j,o1));
        % weights = [weights; matches(i,j)];

        if useMultibody
            % Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o2) / norm(Fs(:,:,i,j,o2));
            % Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,3) / norm(Fs(:,:,i,j,3));
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
[mu_1, sigma_1] = sturm1(Fo, width, height)
[mu_1_1, sigma_1_1] = sturm1_1(Fo, width, height)
[mu_f0, sigma_f0] = sturm2(Fo, width, height)

% fprintf("mu_f0: %f, sigma_f0: %f\n", mu_f0, sigma_f0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = ones(size(Fo,3),1);
weights = weights ./ sum(weights);

Fx = [];
Fy = [];
views = [];

for i = 1:numTrials
    % fprintf("%d\n", i);
    [fx, fy, v] = optimizeFocalLength(Fo, weights, mu_f0, sigma_f0, width, height);
    Fx = [Fx; fx];
    Fy = [Fy; fy];
    views = [views; v];
end

% Remove negative values
% Fx(Fx < 0) = [];
% Fy(Fy < 0) = [];

h = histogram(Fx,min(500, calcnbins(Fx,'fd')));
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);
fx = Fx;
fx(or(fx < th_low, fx > th_high)) = [];
fx = median(fx);

figure;
h = histogram(Fy,min(500, calcnbins(Fy,'fd')));
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);
fy = Fy;
fy(or(fy < th_low, fy > th_high)) = [];
fy = median(fy);

% Find k-NN
I = knnsearch([Fx Fy], [fx, fy], 'K', 3);
v = unique(views(I,:));

% Joint refinement - focal length + principal point
[~,~,u,v] = optimizePrincipalPoints(Fo(:,:,v), fx, fy, width, height);

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
    'StepTolerance',1e-10,...
    'FunctionTolerance',1e-10,...
    'MaxIterations',4,...
    'MaxFunctionEvaluations',1e6,...
    'TolFun', 1e-10,...
    'TolX',1e-10);

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