function E = costFunctionKruppaClassicalWeighted(F, X, weights)
%% COSTFUNCTIONKRUPPACLASSICAL computes Classical Kruppa Cost function to find the Optimal Intrinsic Parameters

%   Input
%       F    - Fundamental Matrix between given two Images
%       X    - Approximate Values of Intrinsics (1 x 5)
%
%   Output
%       E    - Computed Cost

%% Function starts here

% Transform Intrinsics to Matrix Form
K = [X(1) 0 1392/2; 0 X(2) 512/2; 0 0 1];

% Compute the Conic W^-1 (Image of the Absolute Conic)
W_inv = K*K';

% Initialize Cost
E = [];

% Compute the Cost using Classical Kruppa's Equation
for i = 1:size(F,3)
    % 1st term of Kruppa's Equation
    A = F(:,:,i) * W_inv * F(:,:,i)';

    A = A/norm(A,'fro'); % Compute Frobenius Norm

    % 2nd term of Kruppa's Equation
    [~,~,V] = svd(F(:,:,i)'); % Compute SVD of Fundamental Matrix

    V = V(:,end); % Last Colum Vector is the Epipole

    Epi = skewSymmetricMatrix(V);

    B = Epi * W_inv * Epi';

    B = B/norm(B,'fro'); % Compute Frobenius Norm

    % Compute Cost
    E1 = A - B;

    % Append Computed Cost
    E = [E E1(1,1:3) E1(2,2:3)];
end

E = E .* weights';

end
