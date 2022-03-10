function [ePPMs H iter] = proj2euc(PPMs,Kest,opts)
%
% Euclidean upgrade using DIAC constraints.
%
%   PPMs    3x4xN   stacked projective PPMs
%   Kest    3x3(xN) instrinsics estimate(s)
%   opts    options structure. See proj2euc_options
%   ePPMs   3x4xN   stacked euclidean PPMs
%   H       4x4     ePPM = PPM * H
%   iter    algorithm iterations

% Author: Riccardo Gherardi (riccardo.gherardi@univr.it)
% Acknowledgements: based on a previous version by A. Fusiello
% References:
%
%   [1] Heyden, Cipolla. "A linear iterative method...", 2001
%   [2] Pollefeys et al. "Surviving dominant planes...", 2002

% TODO - have weights and constraints on length 4, not 6

    if nargin == 0; ePPMs = proj2euc_options(); return; end

    error(nargchk(1,3,nargin))
    error(nargoutchk(0,3,nargout))
    
    if nargin < 2; Kest = eye(3); end
    if nargin < 3; opts = proj2euc_options(); end
    if size(Kest,3) == 1; Kest = repmat(Kest,[1 1 size(PPMs,3)]); end
    if size(PPMs,3) * sum(opts.constraints) < 10; error('Not enough cameras.'); end
    
    ePPMs = PPMs; % preallocation
    pPPMs = PPMs; % preallocation
    for cnt = 1:size(pPPMs,3) % preconditioning
       
        PPMs(:,:,cnt) = PPMs(:,:,cnt) / norm(PPMs(3,1:3,cnt));
        pPPMs(:,:,cnt) = Kest(:,:,cnt) \ PPMs(:,:,cnt); % inv(K)*P
        
    end

    DIAC = zeros(3,3,size(pPPMs,3)); % DIAC estimates
    d = duplication(4); % M&N duplication matrix

    c = opts.constraints; % enabled constraints on w* = K*K'
    ews = ones(size(pPPMs,3),1); % initial relative equation weights (from [1])
    if ~opts.constraintweights; cws = ones(6,1); % contraint weights (from [2])
    else cws = [5 100 10 10 0.1 0.1]; % fx-fy sk x0 y0 fx-1 fy-1
    end
    
	exitf = 1;
    for iter = 1:opts.maxiter
        
        a = []; %#ok<*AGROW> constraint matrix
        for cnt = 1:size(pPPMs,3) % build constraint matrix

            ew = ews(cnt); % equation weight

            if c(1) % fx-fy
                
                a1 = kron(pPPMs(1,:,cnt), pPPMs(1,:,cnt));
                a1 = a1 - kron(pPPMs(2,:,cnt), pPPMs(2,:,cnt));
                a = [a; a1 * d * cws(1) / ew];
                
            end

            if c(2); a = [a; kron(pPPMs(2,:,cnt), pPPMs(1,:,cnt)) * d * cws(2) / ew]; end % sk
            if c(3); a = [a; kron(pPPMs(3,:,cnt), pPPMs(1,:,cnt)) * d * cws(3) / ew]; end % x0
            if c(4); a = [a; kron(pPPMs(3,:,cnt), pPPMs(2,:,cnt)) * d * cws(4) / ew]; end % y0

            if c(5) % fx-1
                
                a2 = kron(pPPMs(1,:,cnt), pPPMs(1,:,cnt));
                a2 = a2 - kron(pPPMs(3,:,cnt), pPPMs(3,:,cnt));
                a = [a; a2 * d * cws(5) / ew];

            end
            
            if c(6) % fy-1
                
                a3 = kron(pPPMs(2,:,cnt), pPPMs(2,:,cnt));
                a3 = a3 - kron(pPPMs(3,:,cnt), pPPMs(3,:,cnt));
                a = [a; a3 * d * cws(6) / ew];
                
            end

        end
        
        DAQ = ns(a);
        DAQ = ivech(DAQ);
        DAQ = DAQ' + tril(DAQ,-1); % DAQ, omega star
        
        DAQ = (DAQ + DAQ') / 2; % symmetry
        [U D V] = svd(DAQ);
        D(D < 0) = 0; % nearest posdef in frob^2
        D(4,4) = 0; % rank 3 approximation
        DAQ = U * D * V';
        DAQ = DAQ / DAQ(3,3); % normalization
        
        meanf = 0; % average focal lenght
        for cnt = 1:size(pPPMs,3) % compute DIACs
            
            DIAC(:,:,cnt) = pPPMs(:,:,cnt) * DAQ * pPPMs(:,:,cnt)';
            if opts.equationweights; ews(cnt) = DIAC(3,3,cnt); end
            DIAC(:,:,cnt) = DIAC(:,:,cnt) / DIAC(3,3,cnt);

%             fx = DIAC(1,1,cnt) - DIAC(1,3,cnt)^2;
%             fy = DIAC(2,2,cnt) - DIAC(2,3,cnt)^2;
%             
%             meanf = meanf + real(sqrt(fx)); assert(imag(fx) == 0)
%             meanf = meanf + real(sqrt(fy)); assert(imag(fy) == 0)
            
            % TODO - if opts.usemeanfocal
            % TODO - or if equationweights

            fx = abs(DIAC(1,1,cnt) - DIAC(1,3,cnt)^2);
            fy = abs(DIAC(2,2,cnt) - DIAC(2,3,cnt)^2);
            
            meanf = meanf + sqrt(fx);
            meanf = meanf + sqrt(fy);
            
        end
        meanf = meanf / (2 * size(pPPMs,3));
        
        prev = exitf;
        exitf = 0;
        for cnt = 1:size(pPPMs,3) % update PPMs

            if ~opts.usemeanfocal; Kcnt = eye(3);
            else Kcnt = diag([meanf meanf 1]);
            end
            Kcnt(1,3) = DIAC(1,3,cnt);
            Kcnt(1,3) = DIAC(2,3,cnt);
            
            Kest(:,:,cnt) = Kest(:,:,cnt) * Kcnt;
            pPPMs(:,:,cnt) = Kest(:,:,cnt) \ PPMs(:,:,cnt);
            exitf = exitf + norm(eye(3) - Kcnt,'fro');

        end        
        exitf = exitf / size(pPPMs,3);
        if exitf < opts.tol && exitf > prev; break; end
        
        if opts.debug % show current results
                
            [U D] = svd(DAQ);
            H = [U * sqrt(D(1:4,1:3)), [0 0 0 1]'];
            for cnt = 1:size(pPPMs,3) % compute Ks
                
                K = KG(PPMs(:,:,cnt) * H);
                fprintf('%6.2f  ',vech(K'));
                fprintf('\n');
                
            end
            fprintf('\neye - Kest = %f\n\n',exitf);
        
        end

    end

    [U D] = svd(DAQ);
    H = [U * sqrt(D(1:4,1:3)), [0 0 0 1]'];
    for cnt = 1:size(PPMs,3); ePPMs(:,:,cnt) = PPMs(:,:,cnt) * H; end
    if iter == opts.maxiter; warning('autocal:maxiter','Max iterations reached.'); end

function opt = proj2euc_options()

    opt.tol = 1e-2;
    opt.debug = false;
    opt.maxiter = 100;
    opt.usemeanfocal = false;
    opt.equationweights = true;
    opt.constraintweights = true;
    opt.constraints = [1 1 0 0 0 0];
