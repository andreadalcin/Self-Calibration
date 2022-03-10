clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')

% Load camera parameters
% K - intrinsic parameters (gt)
% width - image width
% height - image height
load('Synthetic/Data/camera_canon.mat')

% Parameters
num_cams = 10;
motions = 1;

sigma = 0;
num_samples = 100;
min_outlier_f = 1e1;
max_outlier_f = 1e5;
outlier_ratio = 0.3;

matches = zeros(1,200*motions,4);
index = 1;

num_outliers = nchoosek(num_cams, 2) * outlier_ratio;

o = 0;

for i = 1:num_cams-1
    for j = i+1:num_cams
        if o < floor(num_outliers)
            r = (max_outlier_f - min_outlier_f) .* rand(1) + min_outlier_f;
            K_curr = [r 0 width / 2; 0 r height / 2; 0 0 1];
            o = o + 1;
        else
            K_curr = K;
        end

        matches(index,:,:) = sampleMatches(K_curr, sigma, motions);
        index = index + 1;
    end
end

save(sprintf('/Users/andreadalcin/Documents/GitHub/MBSfM-Autocalibration/matches/mathces-%d-%d.mat', num_cams, motions), 'matches')