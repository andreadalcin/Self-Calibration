function [rPPMs H ePPMs] = selfCal(pPPMs,Vs,opts)
%
% function H = newAutocal(pPPMs,Vs,opts)
%
%   Will return a euclidean upgrade of the perspective projection
%   matrices pPPMs, that is ePPMs(:,:,k) = pPPMs(:,:,k) * H. Vs contains
%   the viewport matrices of each camera, which can be computed from
%   the corresponding image size.
%
%       opts = newAutocal;
%       ePPMs = newAutocal(pPPMs,[1024 768]);
%       [ePPMs H] = newAutocal(pPPMs,Vs,opts);
% 
% Author: Riccardo Gherardi (riccardo.gherardi@univr.it)

    if nargin == 0; ePPMs = options(); return; end
    
    error(nargchk(2,3,nargin))
    error(nargoutchk(1,3,nargout))
    
    if nargin < 3; opts = options(); end
    if numel(Vs) == 2 % width, heigth
        
        w = Vs(1); % width
        h = Vs(2); % height
        Vs = norm([w h]) * eye(3);
        Vs(:,3) = [w h 2] / 2;
        
    end
    if size(Vs,3) == 1; Vs = repmat(Vs,[1 1 size(pPPMs,3)]); end

    % --- choose cost operator
    
    costop = @(x) realpow(x,2);
    if opts.absORsquare; costop = @abs; end

    % --- normalize Ps
    
    T = [Vs(:,:,1) \ pPPMs(:,:,1); 0 0 0 1];
    for cnt = 1:size(pPPMs,3) % ensure pPPMs(:,:,1) == eye(3,4)
        
        pPPMs(:,:,cnt) = Vs(:,:,cnt) \ pPPMs(:,:,cnt) / T; % normalize
        pPPMs(:,:,cnt) = pPPMs(:,:,cnt) / norm(pPPMs(3,1:3,cnt));
    
    end
    
    % --- iterate over focal pairs
    
    focalrange = createRange(opts.focalmethod,opts.samples);
    hwaitbar = waitbar(0,'Searching focals','Name','Self calibration');
    costs = zeros(opts.samples,opts.samples,size(pPPMs,3) - 1,4);
    w = opts.weights;
    
    tic
    prevtime = 0;
    P2 = pPPMs(:,:,2);
    for cnt2 = 1:opts.samples

        f2 = focalrange(cnt2);
        K2 = diag([f2 f2 1]);
        G2 = K2 \ P2;
        t2 = G2(:,4);
        
        r1 = t2 / norm(t2);
        r2 = cross(r1,opts.v);
        r3 = cross(r1,r2);
        r3 = r3 / norm(r3);
        r2 = r2 / norm(r2);
        R = [r1 r2 r3]; % R' * t2 = [norm(t2) 0 0]'        

        for cnt1 = 1:opts.samples

            f1 = focalrange(cnt1);
            K1 = diag([f1 f1 1]);
            R2 = G2(:,1:3) * K1;

            R2 = R' * R2;
            scale = norm(R2(3,:));
            R2 = R2 / scale; % absorbs scale ambiguity
            r1 = cross(R2(2,:),R2(3,:));
            pinf = r1 - R2(1,:);
            pinf = pinf / norm(t2) * scale;

            H = [K1 [0 0 0]'; pinf 1];
            
%             for cnt3 = 1:size(costs,3)
%                 
%                 K = KG(pPPMs(:,:,cnt3 + 1) * H);
%                 
%                 c1 = costop(w(1) * (K(1,1) - K(2,2)));   % fx - fy
%                 c2 = costop(w(2) *  K(1,2));             % skew
%                 c3 = costop(w(3) *  K(1,3));             % x0
%                 c4 = costop(w(3) *  K(2,3));             % y0
%                 
%                 costs(cnt1,cnt2,cnt3,1) = c1;
%                 costs(cnt1,cnt2,cnt3,2) = c2;
%                 costs(cnt1,cnt2,cnt3,3) = c3;
%                 costs(cnt1,cnt2,cnt3,4) = c4;
% 
%                 % after verrifying everything works:
%                 % costs(cnt1,cnt2) = costs(cnt1,cnt2) + c1..4;
%                 
%             end

            cost = 0;
            for cnt3 = 1:size(costs,3)
                
                K = art(pPPMs(:,:,cnt3 + 1) * H); % KG
                
                cost = cost + costop(w(1) * (K(1,1) - K(2,2)));   % fx - fy
                cost = cost + costop(w(2) *  K(1,2));             % skew
                cost = cost + costop(w(3) *  K(1,3));             % x0
                cost = cost + costop(w(3) *  K(2,3));             % y0
                
            end
            costs(cnt1,cnt2,1,1) = cost;

        end
        
        if prevtime ~= floor(toc) % update waitbar once a second
        
            str = sprintf('Elapsed time: %.3g',toc);
            waitbar(cnt2 / opts.samples, hwaitbar, str);
        
        end
        prevtime = floor(toc);
    
    end
    fprintf('Elapsed time: %.3g\n',toc);
    close(hwaitbar)

    % --- print debug info, aggregate costs

    % should be redone, because assume costs is modified in place
    % if opts.debug; bestSolutionForConstraint(Vs,focalrange,costs); end
    % if opts.debug; plotCostFunctions(costs,cnt1,cnt2,opts); end
    
    costs = sum(costs, 4);
    costs = sum(costs, 3);
    
    % --- find the best solution and recover H

    [cnt1 cnt2] = find(costs == min(costs(:)));
    
    f2 = focalrange(cnt2);
    K2 = diag([f2 f2 1]);
    G2 = K2 \ P2;
    t2 = G2(:,4);
        
    r1 = t2 / norm(t2);
    r2 = cross(r1,opts.v);
    r3 = cross(r1,r2);
    r3 = r3 / norm(r3);
    r2 = r2 / norm(r2);
    R = [r1 r2 r3]; % R' * t2 = [norm(t2) 0 0]'        

    f1 = focalrange(cnt1);
    K1 = diag([f1 f1 1])
    R2 = G2(:,1:3) * K1;

    R2 = R' * R2;
    scale = norm(R2(3,:));
    R2 = R2 / scale; % absorbs scale ambiguity
    r1 = cross(R2(2,:),R2(3,:));
    pinf = r1 - R2(1,:);
    pinf = pinf / norm(t2) * scale;
    H = [K1 [0 0 0]'; pinf 1];
    
    % --- compute euclidean upgrade. Rejoice.
    
    ePPMs = pPPMs; % preallocation
    for cnt = 1:size(ePPMs,3)
        
        ePPMs(:,:,cnt) = Vs(:,:,cnt) * pPPMs(:,:,cnt) * H;
        ePPMs(:,:,cnt) = ePPMs(:,:,cnt) / norm(ePPMs(3,1:3,cnt));
    
    end

    % --- refine
    
    rPPMs = pPPMs; % preallocation
    for cnt = 1:size(rPPMs,3)
        
        rPPMs(:,:,cnt) = pPPMs(:,:,cnt) * H; % upgraded but not denormalized
        rPPMs(:,:,cnt) = rPPMs(:,:,cnt) / norm(rPPMs(3,1:3,cnt));
    
    end
    
    x = [1 0 0 0 0 0 1]; % f x0 y0 pinf(1:4)
    opt = optimset('lsqnonlin');
    opt = optimset(opt,'algorithm','levenberg-marquardt');
    if ~opts.debug; opt = optimset(opt,'display','off');
    else opt = optimset(opt,'display','iter-detailed');
    end
    x = lsqnonlin(@(x) costfun(rPPMs,x,w),x,[],[],opt);
    Href = HfromX(x);
    
    for cnt = 1:size(rPPMs,3)
        
        rPPMs(:,:,cnt) = Vs(:,:,cnt) * rPPMs(:,:,cnt) * Href;
        rPPMs(:,:,cnt) = rPPMs(:,:,cnt) / norm(rPPMs(3,1:3,cnt));
    
    end
    
function H = HfromX(vH)

    H = diag([vH(1) vH(1) 1 1]); % fx fy
    H(4,:) = vH(4:7); % plane at infty
    H(1:2,3) = vH(2:3); % x0 y0
    
function c = costfun(P,x,w)

    H = HfromX(x);
    c = []; %#ok<*AGROW>
    for cnt = 1:size(P,3)

        K = art(P(:,:,cnt) * H); % KG
        c = [c w(1) * (K(1,1) - K(2,2))];   % fx-fy
        c = [c w(2) *  K(1,2)];             % sk
        c = [c w(3) *  K(1,3)];             % x0
        c = [c w(3) *  K(2,3)];             % y0
%       c = [c w(4) * (K(1,1) - K(3,3))];   % fx-1
%       c = [c w(4) * (K(2,2) - K(3,3))];   % fy-1
        
    end % for cnt

function range = createRange(method,dim)
%
% Creates a range of focal values from 0.3 to 3 inclusive.
% The methods differ in how the values are distributed; in practice
% however I never noticed any advantage in either one of them. In
% particular, powerlaw tries to allocate more samples to more
% probable configurations.
%
% When generating debug plots, its better to not use powerlaw
% method since the resulting graphs are deformed, and therefore
% less intuitive.

    if nargin < 1; method = 'powerlaw'; end
    if nargin < 2; dim = 50; end
    switch lower(method)
    
        case 'linear', range = linspace(0.3,3,dim);
        case 'logarithmic', range = 3 * logspace(-1,0,dim);
        otherwise % powerlaw
    
            r = mod(dim,2);
            left = (((dim/2):-1:1) / (dim/2)) .^ 2;
            right = ((0:((dim-r)/2-1+r)) / ((dim-r)/2-1+r)) .^ 2;
            range = 1 + [-0.7 * left 2 * right];
            
    end
    
function opts = options()
%
% Different v will result in different R, however that does
% not influences the computation of the plane at infinity.
% Higher number of samples will result in longer execution
% times (possibly higher precision).

    opts.debug = false;                 % enables debugging output
    opts.v = [0 0 1]';                  % direction used in creating R
    opts.samples = 70;                  % f1,f2 samples in [0.3..3]
    opts.focalmethod = 'logarithmic';   % linear, logarithmic, powerlaw
    opts.absORsquare = true;            % true = abs, false = square
    opts.saveImages = false;            % save cost profiles to disk
    % opts.weights = [1 1 0 0];         % standard
    % opts.weights = [100 100 10 1];    % non standard
    % opts.weights = [5 100 0 0];       % half pollefeys
    opts.weights = [5 100 10 0.1];      % pollefeys weights

function plotCostFunctions(costs,row,col,opts)
%
%   Invoked when debugging. Prints nice figures.

    hardlimits = [0.2 0.01 0.1 0.1]; % ar sk x0 y0
    titles = {'ar','sk','x0','y0'};

    sumcosts = sum(costs,3);
    maxcosts = max(costs,[],3);
    costs(:,:,end + 1,:) = sumcosts; % aggregated costs sum
    costs(:,:,end + 1,:) = maxcosts; % aggregated costs max

    rep = 1; % colormap wraps
    for cnt = 1:size(costs,3)

        figure
        set(gcf,'Name',sprintf('K%d',cnt + 1))
        if cnt == size(costs,3) - 1; set(gcf,'Name','Sum'); end
        if cnt == size(costs,3) - 0; set(gcf,'Name','Max'); end
        
        % subplot(2,2,1), imshow(costs(:,:,cnt,1),[]), colorbar, title('fx - fy')
        % subplot(2,2,2), imshow(costs(:,:,cnt,2),[]), colorbar, title('sk')
        % subplot(2,2,3), imshow(costs(:,:,cnt,3),[]), colorbar, title('x0')
        % subplot(2,2,4), imshow(costs(:,:,cnt,4),[]), colorbar, title('y0')
        
        for cnt2 = 1:4
            
            % subplot(1,4,cnt2)
            subplot(2,2,cnt2)
            c = costs(:,:,cnt,cnt2);
            c(c > hardlimits(cnt2)) = hardlimits(cnt2);
            
            % c(row,col + (-3:3)) = 3 / 4 * hardlimits(cnt2); % draw cross
            % c(row + (-3:3),col) = 3 / 4 * hardlimits(cnt2); % yellow
            
            imshow(c,[])
            hold on
            plot(row,col,'m*')
            colorbar
            title(titles{cnt2})
            
            if opts.saveImages
                
                fname = sprintf('figures/costK%d_%s.png',cnt,titles{cnt2});
            
                % imwrite(c,fname)

                h = figure;
                imshow(c,[])
                hold on
                plot(row,col,'m*')
                colormap(jet)
                colorbar
                saveas(gcf,fname,'png')
                close(h)

            end
                        
        end
        colormap(repmat(jet(256),rep,1))
                
    end

function bestSolutionForConstraint(Vs,focalrange,costs)
%
%   Invoked when debugging.
%   Computes the best solution for each K entry and their sum.
        
    fprintf('\n')
    for cnt = 1:4

        cost = costs(:,:,end,cnt);
        [r c] = find(cost == min(cost(:)));
        
        r = Vs(1,1,1) * focalrange(r);
        c = Vs(1,1,2) * focalrange(c);
        fprintf('Best solution (%d) : %g %g\n',cnt,r,c);
    
    end
    
    cost = sum(costs(:,:,end,:),4); % are they commensurable?
    [r c] = find(cost == min(cost(:))); % are we adding apples and pears?

    r = Vs(1,1,1) * focalrange(r);
    c = Vs(1,1,2) * focalrange(c);
    fprintf('Best solution sum : %g %g\n',r,c);
