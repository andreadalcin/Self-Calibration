function [thisEstimate] = robustlsqnonlin(Fs, ds, X0, Options)

weightMethod = 'cauchy';
convergenceThreshold = 1e-6;

[weightFunction, tuningConstant] = weightFunAndConstant(weightMethod);

thisEstimate = X0;
hasConverged     = false;
previousEstimate = inf(size(X0));
weights          = ones(size(Fs,3),1);
iterationCounter = 1;
while ~hasConverged && iterationCounter < 1e2
    %%% weighted LSQ
    weightedFun = @(X) costFunctionMendoncaCipollaWeighted(Fs, ds, X, weights, '2');

    thisEstimate = lsqnonlin(weightedFun, X0, [], [], Options)

    % thisEstimate = varargout{1};
    hasConverged = norm(thisEstimate - previousEstimate)^2 < convergenceThreshold;
    
    %%% update weights
    residuals = costFunctionMendoncaCipolla(Fs, ds, thisEstimate, '2');
    residuals = residuals(:);
    
    residualLeverages = leverage(residuals);
    robustVar         = mad(residuals, 1);
    
    r = residuals ./ (tuningConstant * robustVar * sqrt(1 - residualLeverages));
    
    % weights = weightFunction(r);
    previousEstimate = thisEstimate;
    iterationCounter = iterationCounter + 1;
end
end

function [weightFun, tuningConstant] = weightFunAndConstant(method)
switch lower(method)
    case 'bisquare'
        weightFun      = @(r) (abs(r) < 1) .* (1 - r.^2).^2;
        tuningConstant = 4.685;
        
    case 'andrews'
        weightFun      = @(r) (abs(r) < pi) .* sin(r) ./ r;
        tuningConstant = 1.339;
        
    case 'cauchy'
        weightFun      = @(r) 1 ./ (1 + r.^2);
        tuningConstant = 2.385;
        
    case 'fair'
        weightFun      = @(r) 1 ./ (1 + abs(r));
        tuningConstant = 1.4;
        
    case 'huber'
        weightFun      = @(r) 1 ./ max(1, abs(r));
        tuningConstant = 1.345;
        
    case 'logistic'
        weightFun      = @(r) tanh(r) ./ r;
        tuningConstant = 1.205;
        
    case 'ols'
        weightFun      = @(r) ones(size(r));
        tuningConstant = 1;
        
    case 'talwar'
        weightFun      = @(r) 1 * (abs(r) < 1);
        tuningConstant = 2.795;
        
    case 'welsch'
        weightFun      = @(r) exp(-(r.^2));
        tuningConstant = 2.985;
        
end
end