function [P1 P2] = PfromF(F,m1,m2,v)

% We expect a F computed from normalized coordinates (wrt image plane).
% This way, we can approximate K1, K2 as identities and Q2 as a rotation.
% The resulting couple should be quasi-euclidean after cheirality.
% In any case, this is at least better than the canonical couple, which
% is guaranteed to have the second camera center at infinity.

% Author: Riccardo Gherardi (riccardo.gherardi@univr.it)

    error(nargchk(1,4,nargin))
    error(nargoutchk(0,3,nargout))
    
    F = F / norm(F);
    [U D V] = svd(F');
    e1 = V(:,3); % left epipole, already normal

    if nargin < 4; v = [0 0 1]'; end % any vector ~= e1 will do
    r2 = cross(e1,v); % different v will result in different Ps
    r3 = cross(e1,r2);
    r3 = r3 / norm(r3);
    r2 = r2 / norm(r2);
    R = [e1 r2 r3]; % R * [1 0 0]' = e1
    
    P1 = eye(3,4); % canonical cameras
    P2 = [star(e1) * F e1];
    
    num = 1; % not needed imho
    for cnt = 1:num
    
        Q2 = P2(:,1:3);
        Q2 = R' * Q2;
        Q2 = Q2 / norm(Q2(3,:));
        Q2(1,:) = cross(Q2(2,:),Q2(3,:));
        Q2(1,:) = Q2(1,:) / norm(Q2(1,:)); % * norm(Q2(2,:));
        Q2 = R * Q2;
        
        delta = 1; % should affect only scale
        P2 = [Q2 delta * e1];

        P = P1;
        P(:,:,2) = P2;
        x = [1 0 0 0 0 0 1]; % f x0 y0 pinf(1:4)
        opt = optimset('lsqnonlin');
        opt = optimset(opt,'algorithm','levenberg-marquardt');
        opt = optimset(opt,'display','off'); % 'iter-detailed');
        x = lsqnonlin(@(x) costfun(P,x),x,[],[],opt);
        H = HfromX(x);
        P1 = P1 * H;
        P2 = P2 * H;
        
        % P1, P2, KG(P2), ...
        % svd(KG(P2)' * F * KG(P1))', ...
        % max(max(abs(F) - abs(fund(P1,P2)))) %#ok<*NOPRT>
               
    end % for cnt = 1:num
    
    if nargin > 1 % cheirality

        if nargin == 2

            m2 = m1(:,:,2);
            m1 = m1(:,:,1);

        end

        K2 = KG(P2);
        [K1 G1] = KG(P1);
        
        [R t fail] = sr(F,K1,K2,m1,m2);
        P2 = K2 * [R t] * [G1; 0 0 0 1];

        if fail; error('Cheirality failed!'); end

    end % if cheirality

    if nargout < 2; P1(:,:,2) = P2; end
    
function H = HfromX(vH)

    H = diag([vH(1) vH(1) 1 1]); % fx fy
    H(4,:) = vH(4:7); % plane at infty
    H(1:2,3) = vH(2:3); % x0 y0
    
function c = costfun(P,x)

    w = [1 1 0 0]; % standard
    % w = [100 100 10 1]; % non standard
    % w = [5 100 0 0]; % half pollefeys
    % w = [5 100 10 0.1]; % pollefeys weights

    H = HfromX(x);
    c = []; %#ok<*AGROW>
    for cnt = 1:size(P,3)

        K = abs(KG(P(:,:,cnt) * H));
        c = [c w(1) * (K(1,1) - K(2,2))];   % fx-fy
        c = [c w(2) *  K(1,2)];             % sk
        c = [c w(3) *  K(1,3)];             % x0
        c = [c w(3) *  K(2,3)];             % y0
        c = [c w(4) * (K(1,1) - K(3,3))];   % fx-1
        c = [c w(4) * (K(2,2) - K(3,3))];   % fy-1
        
    end % for cnt

% should be a tautology in the cost function
%
%     P1 = P(:,:,1) * H;
%     P2 = P(:,:,2) * H;
%     K1 = KG(P1);
%     K2 = KG(P2);
%     F = fund(P1,P2);
%     E = K2' * F * K1;
%     D = svd(E);
% 
%     c = [c D(1) - D(2)]; % D(3)];

% old wrong broken version
%
%         Q2 = P2(:,1:3);
%         Q2 = R' * Q2;
%         Q2(1,:) = cross(Q2(2,:),Q2(3,:));
%         r2 = cross(Q2(3,:),Q2(1,:));
%         fest = Q2(2,:) * r2' / norm(r2); % ||Q2(2,:)|| projected on r2
%         Q2(1,:) = Q2(1,:) / norm(Q2(1,:)) * fest;
%         Q2 = R * Q2;
%         Q2 = Q2 / norm(Q2(3,:));

% some other broken code
%
%         m1 = K1 \ [m1; ones(1,size(m1,2))];
%         m2 = K2 \ [m2; ones(1,size(m2,2))];
%         m1 = [m1(1,:) ./ m1(3,:); m1(2,:) ./ m1(3,:)];
%         m2 = [m2(1,:) ./ m2(3,:); m2(2,:) ./ m2(3,:)];
