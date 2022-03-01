clc, clear, close all
format longG;

% rng default;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

% Generate ground truth parameters
width = 4000;
height = 3000;
f_gt = 2 * (width + height);
K = [f_gt 0 width / 2; 0 f_gt height / 2; 0 0 1];

% num_F = 15;
% num_F = 25;
% num_F = 50;
num_F = 75;

% outlier_ratio = 0;
xaxis_limit = true;

% Sample fundamental matrices
Fs = [];
sigma = 5;
outlier_indices = [];

for i = 1:num_F
    if (rand > 0)
        K_curr = K;
    else
        a = 1e1;
        b = 1e5;
        r = (b-a) .* rand(1) + a;
        K_curr = [r 0 width / 2; 0 r height / 2; 0 0 1];
        outlier_indices = [outlier_indices; i];
    end

    Fs(:,:,i) = sampleFundamental(K_curr,sigma);
end

init_focal_length(Fs, width, height, f_gt, xaxis_limit, num_F);


function [f0] = init_focal_length(Fs, width, height, f_gt, xaxis_limit, num_f)

% Output from methods
[mu_1, sigma_1] = sturm1(Fs, width, height)
[mu_1_1, sigma_1_1] = sturm1_1(Fs, width, height)
[mu_f0, sigma_f0, x, ySix] = sturm2(Fs, width, height);
mu_f0
sigma_f0

% Plot KDE
figure('Name', 'Sturm plot results', 'NumberTitle', 'off');
title(sprintf('Robust initialization - Synthetic test - %d fundamental matrices', num_f))
hold on
h = plot((x - f_gt) / f_gt, ySix ./ max(ySix),'k-','LineWidth',2);
set(h,'LineSmoothing','On')
xlabel('Relative Error (%)');
ylabel('Density');
if xaxis_limit
    xlim([-1 1])
end
xline((mu_1_1 - f_gt) / f_gt,'r-','LineWidth',3);
xline((mu_1 - f_gt) / f_gt,'g-','LineWidth',3);
xline((mu_f0 - f_gt) / f_gt,'b-','LineWidth',3);
hold off

end