function [f, uv] = runTest(X0, Fs, K_approx)

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionKruppaSimplified2(Fs, X),X0,[],[],Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];
disp('Intrinsics computed by Simplified Kruppa`s Self Calibration')
disp(K_SK)

% Print error metrics
disp("focal length")
f = abs(K_SK(1,1) - K_approx(1,1)) + abs(K_SK(2,2) - K_approx(2,2));
disp(f)

disp("principal point")
uv_SK = [K_SK(1,3), K_SK(2,3)];
uv_approx = [K_approx(1,3), K_approx(2,3)];
uv = norm(uv_SK - uv_approx);
disp(uv)

end