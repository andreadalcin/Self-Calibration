function [ Q, K ] = refitHypothesis( X, P, H, epsilon, model, voting)

%REFIT_PREF refit preferences, accept a refit only if the CS is increased
%   at the moment work only for hard preferences
%   INPUT:
%          X            data points
%          epsilon      inlier threshold
%          P            preference matrix
%          H            sampled hypotheses
%          voting       voting function as defined in prefMat 0 hard
%   OUTPUT:
%          Q            new preference matrix
%          K            new sampled hypotheses

if nargin <6
    voting =0;
    disp('refit_consensus: hard voting employed');
end
[ distFun, hpFun,fit_model,  cardmss, ~, ~] = set_model( model );

m   = size(P,2); % number of hypotheses

count = 0;
if (numel(epsilon)==1)
    epsilon=epsilon*ones(1,m);
end
Q = P;
K = H;
for j = 1:m
    % compute consensus set
    inlier  = P(:,j)>0;
    card_cs = sum(inlier); % cardinality of the consensus set
    % refit model
    if(card_cs>cardmss)
        h_new  = fit_model(X(:, inlier)); % new hypothesis
        % update if necessary
        r = res_line(X, h_new);
        if(sum(r<epsilon(j)) >= card_cs)
            K(:,j) = h_new;
            Q(:,j) = prefMat(r, epsilon(j), voting );
            count  = count+1;
        end
        
    end
end

%fprintf('Number of updated hypotheses %i\n', count);
end

