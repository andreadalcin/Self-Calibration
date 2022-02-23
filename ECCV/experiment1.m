clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

load('Dataset/Fs.mat')
load('Dataset/params.mat')

% Fx = [num_cameras, SB, MB]
% Fy = [num_cameras, SB, MB]

sample_size = 1;
num_trials = 100;

width = 4272;
height = 2848;

Fx = [];
Fy = [];

for num_cameras = 3:3
    Fx_sb = [];
    Fy_sb = [];
    Fx_mb = [];
    Fy_mb = [];

    for i = 1:sample_size
        cameras = sort(datasample(1:size(Fs,3), num_cameras, 'Replace', false))

        [fx, fy] = selfCalibrate(false, cameras, Fs, ...,
            pointMatchesInliers1, pointMatchesInliers2, ...
            width, height, num_trials);
        Fx_sb = [Fx_sb; fx];
        Fy_sb = [Fy_sb; fy];
        close all

        [fx, fy] = selfCalibrate(true, cameras, Fs, ...
            pointMatchesInliers1, pointMatchesInliers2, ...
            width, height, num_trials);
        Fx_mb = [Fx_mb; fx];
        Fy_mb = [Fy_mb; fy];
        close all
    end

    disp("fx - SB")
    disp(median(Fx_sb))
    disp("fx - MB")
    disp(median(Fx_mb))
    disp("fy - SB")
    disp(median(Fy_sb))
    disp("fy - MB")
    disp(median(Fy_mb))

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


3,%%
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





function [fx, fy] = selfCalibrate(useMultibody, cameras, Fs, pointMatchesInliers1, pointMatchesInliers2, width, height, numTrials)

num_cameras = size(cameras,2);
p = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load fundamental matrices and priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fo = [];
priors = [];

for a = 1:num_cameras-1
    for b = a+1:num_cameras
        i = cameras(a);
        j = cameras(b);

        if rand > p
            o1 = 1;
            o2 = 2;
            matches1 = pointMatchesInliers1;
            matches2 = pointMatchesInliers2;
        else
            o1 = 2;
            o2 = 1;
            matches1 = pointMatchesInliers2;
            matches2 = pointMatchesInliers1;
        end

        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o1) / norm(Fs(:,:,i,j,o1));
        priors = [priors; matches1(i,j)];
        if useMultibody
            Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o2) / norm(Fs(:,:,i,j,o2));
            priors = [priors; matches2(i,j)];
            Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,3) / norm(Fs(:,:,i,j,3));
            % priors = [priors; pointMatchesInliers3(i,j)];
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
f = [];
for i = 1:size(Fo,3)
    f0 = getInitialEstimateOfFocalLength(Fo(:,:,i), width, height);
    f = [f; f0];
end

% % Filter focal lengths not in the peak bin
h = histogram(f, 50);
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);

f_pct = f;
f_pct(or(f < th_low, f > th_high)) = [];

f0 = median(f_pct);

% Kernel voting
[x,xi] = ksdensity(f,'Bandwidth',f0*0.01,'Kernel','epanechnikov');
[~,I] = max(x);
f0 = xi(I);

fprintf("f0: %f\n", f0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = ones(size(Fo,3),1);
weights = weights ./ sum(weights);

Fx = [];
Fy = [];

for i = 1:numTrials
    % fprintf("%d\n", i);
    [fx_sample, fy_sample] = optimizeParams(Fo, weights, f0, width, height);
    Fx = [Fx; fx_sample];
    Fy = [Fy; fy_sample];
end

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

end


function [fx,fy] = optimizeParams(Fs, weights, f0, width, height)

K0 = [f0, 0, width / 2; 0, f0, height / 2; 0, 0, 1];

Options = optimoptions('lsqnonlin','Display','off', ...
    'StepTolerance',1e-20,...
    'FunctionTolerance',1e-20,...
    'MaxIterations',1e2,...
    'MaxFunctionEvaluations',1e6,...
    'TolFun', 1e-20,...
    'TolX',1e-20);

X0 = [K0(1,1) K0(2,2)];

sample = datasample(1:size(Fs,3),3,'Replace',false,'Weights',weights);
subset = Fs(:,:,sample);

% Randomize principal point
u0 = normrnd(width / 2, width / 6);
v0 = normrnd(height / 2, height / 6);

loss = @(X) costFunctionMendoncaCipollaFocalOnly(subset, X, u0, v0, '2');

K_SK = lsqnonlin(loss, X0, [], [], Options);
K_SK = [K_SK(1) 0 width/2; 0 K_SK(2) height/2; 0 0 1];

fx = K_SK(1,1);
fy = K_SK(2,2);

end