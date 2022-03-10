function [ok] = isMergeableSL(P, L, i, j)
% Check if two clusters A and B can be merged according to single linkage
% The test performs the following steps:
% i) the closest pair (a,b) a in A, b in B is determined
% iia) the jaccard ball Ba = B(a,radius) centered in a is computed
% iib) the jaccard ball Bb = B(b,radius) centered in b is computed
% iii) the intersection Pa and Pb of preferences in Ba and in Bb is computed
% A cluster is mergeable if Pa and Pb have some common element
% If clusters can be merged ok = true; otherwise ok = false;
% 
%  Input: 
%   - P: the preference matrix
%   - L: vector of clusters labels
%   - i and j: index of points belonging to the cluster to be merged
%   - radius: the radius of the jaccard ball
%   - idx: auxiliary variable of points idx
%  Output:
%   - ok: tells if the cluster can be merged 


isInCi = L==L(i);
isInCj = L==L(j);
Pa = min(P(isInCi,:),[],1);
Pb = min(P(isInCj,:),[],1);

% compute jaccard index
if(all(Pa==0))
    jacc = 0;
elseif(all(Pb==0))
    jacc = 0;
else
    s = sum(Pa.*Pb);
    jacc = s/(Pa*Pa'+ Pb*Pb' - s);
end


ok = jacc> 0.0;

end