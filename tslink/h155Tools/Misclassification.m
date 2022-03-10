%--------------------------------------------------------------------------
% This function takes the groups resulted from spectral clutsering and the
% ground truth to compute the misclassification rate.
% groups: [grp1,grp2,grp3] for three different forms of Spectral Clustering
% s: ground truth vector
% Missrate: 3x1 vector with misclassification rates of three forms of
% spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2010
% Modified by Filippo Leveni & Luca Magri to deal with cllusterings with
% different cardinalities

%--------------------------------------------------------------------------


function Missrate = Misclassification(groups,s)
assert(size(groups,2)==1,'groups must be a column vector!')

n = max(max(s),max(groups));
if n > 8
    Missrate(1,1) = 1;
else
    for i = 1:1
        Missrate(i,1) = missclassGroups( groups(:,i),s,n ) ./ length(s);
    end
end
end