clc, clear, close all

load('../K_Canon.mat')
% K = A;

yaw = rand(1) * 30;
roll = rand(1) * 30;
pitch = rand(1) * 30;
R = rotYPR(yaw, roll, pitch);

P1 = [K, [10, 5, 20]'];
P2 = K * [R, [1, 1, 1]'];
F = fund(P1, P2);
% F = sampleFundamental(K, 0);

tau = 1;
u0 = 2136;
v0 = 2136;
f0 = 2e4;

% Undo intrinsic parameters
G = [tau, 0, 0; 0, 1, 0; u0, v0, 1] * F * [tau, 0, u0; 0, 1, v0; 0, 0, 1];
G_1 = diag([f0,f0,1]) * G * diag([f0,f0,1]);
G_1 = G_1 / norm(G_1, 'fro');

% SVD
[U,S,V] = svd(G_1);
u13 = U(1,3);
u23 = U(2,3);
v13 = V(1,3);
v23 = V(2,3);
a = S(1,1);
b = S(2,2);

% Quadratic equation
c1 = a^2 * (1 - u13^2) * (1 - v13^2) - b^2 * (1 - u23^2) * (1 - v23^2);
c2 = a^2 * (u13^2 + v13^2 - 2 * u13^2 * v13^2) - b^2 * (u23^2 + v23^2 - 2 * u23^2);
c3 = a^2 * u13^2 * v13^2 - b^2 * u23^2 * v23^2;

x = zeros(2,1);
d = sqrt(c2^2 - 4*c1*c3);
x(1) = ( -c2 + c3 ) / (2*c1);
x(2) = ( -c2 - c3 ) / (2*c1);

f = sqrt(x) * f0;
disp("Diff abs(f - f0)")
disp(abs(f - K(1,1)))

disp("Recovered f")
[~,I] = min(abs(f - f0));
disp(f(I))

% Linear equations
f1 = - (u23 + v13 * (a * u13 * v13 + b * u23 * v23)) / (a * u13 * u23 * (1 - v13^2) + b * v13 * v23 * (1 - u23^2));
% sqrt(f1) * f0
f2 = - (u13 + v23 * (a * u13 * v13 + b * u23 * v23)) / (a * v13 * v23 * (1 - u13^2) + b * u13 * u23 * (1 - v23^2));
% sqrt(f2) * f0