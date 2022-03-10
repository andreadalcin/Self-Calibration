function [K G V mul] = KG(PPM,dims,height)

% [K G V mul] = KG(PPM,dims,height)
%
%   Factor a projection matrix into extrinsic G = [R t] and
%   intrinsic K parameters. Focal length is assumed positive.
%
%   PPM = rand(3,4);
%   [K G] = KG(PPM); % basic usage
%   [K G V mul] = KG(PPM,size(img));
%   PPM, m*V*K*G

%   Author: Riccardo Gherardi (riccardo.gherardi@univr.it)

    error(nargchk(1,3,nargin));
    error(nargoutchk(0,4,nargout));
    if ~isequal(size(PPM),[3 4]); error('invalid size'); end
    if rank(PPM)<3; warning('cvlab:rankdeficient','rank deficient PPM'); end

    V = eye(3); % viewport matrix
    if nargin>1
        
        if nargin>2; dims(2) = height; end
        V = eye(3) * norm(dims);
        V(:,3) = [dims / 2 1];
        PPM = V \ PPM; % better conditioned?
        
    end
    
    mul = norm(PPM(3,1:3));
    PPM = PPM / mul;
    
    r3 = PPM(3,1:3);
    r2 = PPM(2,1:3);
    r1 = cross(r2,r3);
    r2 = cross(r3,r1);
    r1 = r1 / norm(r1);
    r2 = r2 / norm(r2);
    R = [r1; r2; r3];
    
    K = triu(PPM(:,1:3) * R'); % triu to fix for near eps values
    G = K \ PPM;
    fix = eye(3) .* sign(K + eps); % eps for cameras with center
    mul = det(fix) * mul;          % on the plane at infinity
    K = K * fix;
    G = det(fix) * fix * G;
    if nargout<3; K = V*K; end
