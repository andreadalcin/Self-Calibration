clear, clc, close all;

format long;

% Fix RNG
rng('default')

% Load the Given Data
load('data.mat');

K_approx = A;
disp('Intrinsics')
disp(K_approx)

% Params
num_F = 20;
sigma_F = 0;
sigma_X0 = 1000;

% Fs
Fs = zeros(3,3,num_F);
for i = 1:size(Fs,3)
    R = sampleR();
    t = sampleT();

    T_skew = [0 -t(3) t(2) ; t(3) 0 -t(1) ; -t(2) t(1) 0 ];

    F = inv(K_approx)' * T_skew * R * inv(K_approx);

    % Noise -- multiply by small rotation
    F_noisy = F + sigma_F * randn(3,3);
    Fs(:,:,i) = F_noisy;
end

% Perturbate initial guess - X0
X0 = [A(1,:) A(2,2:3)];
X0 = X0 + sigma_X0 * randn(1,5);
disp("Initial guess")
disp(X0)

%% Run tests

f_array = zeros(size(Fs,3), 1);
uv_array = zeros(size(Fs,3), 1);

for i = 1:size(Fs,3)
    [f, uv] = runTest(X0, Fs(:,:,1:i), K_approx);

    f_array(i) = f;
    uv_array(i) = uv;
end

figure
plot(f_array)
title("Focal length - error vs number of F")

figure
plot(uv_array)
title("Principal point - error vs number of F")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = sampleR()
R = orth(randn(3,3));
end

function t = sampleT()
t = randn(3,1);
t(3) = 1;
end