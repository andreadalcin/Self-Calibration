clc, clear, close all
format long;
addpath('./Synthetic/')
run("ComputerVisionToolkit/cvt_setup.m");

rng('default')

% Load data
load('Data/K_Canon.mat')
load('Data/K_approx.mat')

sigma = 0;
num_cameras = 10;

FN = zeros(3*num_cameras,3*num_cameras);
pointMatchesInliers = zeros(num_cameras,num_cameras);

for i = 1:num_cameras-1
    for j = i+1:num_cameras
        [F, ~] = sampleFundamental(K, sigma);
        
        FN(3*i-2:3*i,3*j-2:3*j) = F;
        FN(3*j-2:3*j,3*i-2:3*i) = F';

        pointMatchesInliers(i,j) = round(10 * rand(1,1));
    end
end

% Step 2
[finalTriplets] = buildTripletsAviodCollinearFast(pointMatchesInliers,FN);
[ FN ] = normalizeForbineusNorm( FN );
[Xs,Y] = optimizeTripletsIRLS( FN,finalTriplets(1:end,1:3),false,false ); 

