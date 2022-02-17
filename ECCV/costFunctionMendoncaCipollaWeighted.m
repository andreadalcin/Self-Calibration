function E = costFunctionMendoncaCipollaWeighted(F, X, weights, Method)
%% COSTFUNCTIONMENDONCACIPOLLA computes Mendonca & Cipolla Cost function to find the Optimal Intrinsic Parameters

%   Input
%       F      - Fundamental Matrix between given two Images
%       X      - Approximate Values of Intrinsics
%       Method - Cost Function Formula {'1' or '2'}
%
%   Output
%       E    - Computed Cost

%% Function starts here

% Transform Intrinsics to Matrix Form
K = [X(1) X(2) X(3); 0 X(4) X(5); 0 0 1];

% Initialize Cost
E = [];

% For the Denominator term of Mendonca & Cipolla's Equation
% N = size(F,3);

% Den = N*(N-1)/2; % For N Images we can find N(N-1)/2 Fundamental Matrix
Den = size(F,3);

% Compute the Cost using Mendonca & Cipolla's Equation
for i = 1:size(F,3)

    % Compute the Essential Matrix 'EM'
    EM = K' * F(:,:,i) * K;

    % Compute SVD of Essential Matrix
    [~,D,~] = svd(EM);

    % Singular Values (3rd value, D(3,3) is 0)
    r = D(1,1);
    s = D(2,2);

    % Compute Cost
    switch Method
        case '1'
            % Use the  Mendonca & Cipolla's Equation-1
            E1 = 1/Den * (r - s)/s;

        case '2'
            % Use the  Mendonca & Cipolla's Equation-2 (constant intrinsics)
            E1 = 1/Den * (r - s)/(r + s);
            % E1 = 0.5 * 1/Den * (r - s)/(r + s) + 0.5 * (K(1,1) - K(2,2)) / (K(1,1) + K(2,2));
        case '3'
            E1 = 1/Den * (r - s)/(r + s);
    end

    % Append Computed Cost
    E = [E E1];
end

E = E .* weights';

[~,~,~, outliers] = robustcov(E);
if size(outliers,1) == size(E,2)
    E(outliers) = E(outliers) * 0;
else
    disp(size(E))
    disp(size(outliers))
end

end
