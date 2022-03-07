clc, clear, close all
format longG;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

load('Dataset/Sturm/fountain.mat')
load('Dataset/params.mat')

% Load parameters
fx_gt = cameraParams.Intrinsics.FocalLength(1);
fy_gt = cameraParams.Intrinsics.FocalLength(2);
u_gt = cameraParams.Intrinsics.PrincipalPoint(1);
v_gt = cameraParams.Intrinsics.PrincipalPoint(2);
width = cameraParams.ImageSize(2);
height = cameraParams.ImageSize(1);
matches = [];

kNN_support = 10;
num_cameras = size(Fs,3);

% Experiments
sample_size = 1;
num_trials = 2000;

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
            width, height, num_trials, kNN_support);
        Fx_sb = [Fx_sb; fx];
        Fy_sb = [Fy_sb; fy];
        U_sb = [U_sb; u];
        V_sb = [V_sb; v];
        close all

        [fx, fy, u, v] = selfCalibrate(true, cameras, Fs, ...
            width, height, num_trials, kNN_support);
        Fx_mb = [Fx_mb; fx];
        Fy_mb = [Fy_mb; fy];
        U_mb = [U_mb; u];
        V_mb = [V_mb; v];
        close all
    end

    % Single body
    fx = median(Fx_sb)
    fy = median(Fy_sb)
    u = median(U_sb)
    v = median(V_sb)
    err_f_sb = 0.5 * 100 * (abs((fx - fx_gt)/fx_gt) + abs((fy - fy_gt)/fy_gt));
    err_pp_sb = 0.5 * 100 * (abs((u - u_gt)/u_gt) + abs((v - v_gt)/v_gt));

    % Multi body
    fx = median(Fx_mb)
    fy = median(Fy_mb)
    u = median(U_mb)
    v = median(V_mb)
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


function [fx, fy, u, v] = selfCalibrate(useMultibody, cameras, Fs, width, height, numTrials, knnSupport)

num_cameras = size(cameras,2);
p = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fundamental matrices and priors
Fo = [];
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

        if Fs(:,:,i,j,o1) ~= zeros(3,3)
            Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o1) / norm(Fs(:,:,i,j,o1),2);
        end
        
        if useMultibody
            if Fs(:,:,i,j,o2) ~= zeros(3,3)
                Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o2) / norm(Fs(:,:,i,j,o2),2);
            end
%             if Fs(:,:,i,j,3) ~= zeros(3,3)
%                 Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,3) / norm(Fs(:,:,i,j,3));
%             end
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
ED = zeros(size(Fo,3),1);
W = ones(size(Fo,3),1);

iter = 0;
K_prev = [mu_f0, 0, width / 2; 0, mu_f0, height / 2; 0, 0, 1];

while iter < 100
    K_curr = optimize_ls(Fo, W, K_prev);
    dK = K_curr - K_prev;
    if dK == zeros(3,3)
        break;
    end
    K_prev = K_curr;

    for i = 1:size(Fo,3)
        dE = compute_dE(dK, Fo(:,:,i), K_curr);
        ED(i) = compute_ED(dE, Fo(:,:,i), K_curr);
    end

    sj = robstd(ED);
    W = exp(-ED ./ sj);

    iter = iter + 1;
end

fx = K_curr(1,1);
fy = K_curr(2,2);
u = K_curr(1,3);
v = K_curr(2,3);

end


function f = kernel_voting(f, bandwidth_factor)
f(f < 100 | f > 1e5) = [];
bandwidth = median(f) * bandwidth_factor;
pdSix = fitdist(f,'Kernel','Width',bandwidth);
x = min(f):.1:max(f);
ySix = pdf(pdSix,x);
[~,I] = max(ySix);
f = x(I);
end


function K_SK = optimize_ls(Fs, W, K_curr)

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
X0 = [K_curr(1,1) K_curr(2,2), K_curr(1,3), K_curr(2,3)];

loss = @(X) costFunctionMendoncaCipollaWeighted(Fs, X, W, '2');

K_SK = lsqnonlin(loss, X0, [], [], Options);
K_SK = [K_SK(1) 0 K_SK(3); 0 K_SK(2) K_SK(4); 0 0 1];
end


function dE = compute_dE(dK, F, K)
dE = dK'*F*K + K*F*dK' + dK'*F*dK;
end

function ED = compute_ED(dE, F, K)
EM = K' * F * K;
[~,D,~] = svd(EM);
ED = ((D(1,1) - D(2,2)) / (2 * norm(dE, 2)));
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