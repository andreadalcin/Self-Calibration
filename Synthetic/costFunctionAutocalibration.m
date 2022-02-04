clc, clear, close all

load('../Data.mat')
K = A;
K(1,1) = 8e2;
K(2,2) = 8e2;
K(1,3) = 1;
K(2,3) = 1;

% F = sampleFundamental(A, 0);
% F = Fs(:,:,1,4)

P1 = [K, [1, 0, 0]'];
P2 = K * [rotYPR(90, 0, 0), [1, 1, 0]'];
F = fund(P1, P2);

u0 = K(1,3);
v0 = K(2,3);

[U,D,V] = svd(F);

u1 = U(:,1);
u2 = U(:,2);
u3 = U(:,3);

v1 = V(:,1);
v2 = V(:,2);
v3 = V(:,3);

r = D(1,1);
s = D(2,2);

syms f
C = [f + u0^2,  u0 * v0,   u0;
     u0 * v0,   f + v0^2,  v0;
     u0,        v0,        1];

eqn1 = (r^2 * v1' * C * v1) * (-u2' * C * u1) - (u2' * C * u2) * (r * s * v1' * C * v2) == 0;
eqn2 = (r^2 * v1' * C * v1) * (u1' * C * u1) - (u2' * C * u2) * (s^2 * v2' * C * v2) == 0;

% eqn1 = (v2' * C * v2) / (r^2 * u1' * C * u1) - (-v2' * 2 * C * v1) / (s * r * u1' * C * u1) == 0;
% eqn2 = (v2' * C * v2) / (r^2 * u1' * C * u1) - (v1' * C * v1) / (s^2 * u2' * C * u2);

[X] = solve([eqn1], [f]);

fplot([lhs(eqn1) rhs(eqn1)], [-5 5])

double(X)
% double(Y)

sqrt(double(X))
% sqrt(double(Y))