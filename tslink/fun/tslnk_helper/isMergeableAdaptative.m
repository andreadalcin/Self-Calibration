function [ok,Ba,Bb,a,b] = isMergeableAdaptative(Jdist, R, NI, P, L, i, j, idx)
% Check if two clusters A and B can be merged.
% The test performs the following steps:
% i) the closest pair (a,b) a in A, b in B is determined
% iia) the jaccard ball Ba = B(a,radius) centered in a is computed
% iib) the jaccard ball Bb = B(b,radius) centered in b is computed
% iii) the intersection Pa and Pb of preferences in Ba and in Bb is computed
% A cluster is mergeable if Pa and Pb have some common element
% If clusters can be merged ok = true; otherwise ok = false;
%
%  Input:
%   - Jdist: the jaccard distance matrix between points
%   - R: rank of indices per rows sort on distances
%   - NI: nearest index matrix:  NI(u,v) = h iff v is the h-nn of v
%   - P: the preference matrix
%   - L: vector of clusters labels
%   - i and j: index of points belonging to the cluster to be merged
%   - radius: the radius of the jaccard ball
%   - idx: auxiliary variable of points idx
%  Output:
%   - ok: tells if the cluster can be merged
if(nargin<7)
    idx = 1:size(Jdist,1);
end

isInCi = L==L(i);
isInCj = L==L(j);
% find closest pair (a,b)
iIsSingleton = sum(isInCi)==1;
jIsSingleton = sum(isInCj)==1;
if(~iIsSingleton || ~jIsSingleton)
    [a,b] = findClosestInClusters(Jdist, isInCi,isInCj);
else
    a = i;
    b = j;
end
assert(L(a)==L(i));
assert(L(b)==L(j))

% find neighbourhood adaptatively
ha = NI(a,b);
hb = NI(b,a);

Ba = R(a, L(R(a,:))==L(a));
if(numel(Ba)>ha)
    Ba = Ba(1:ha);
end

Bb = R(b, L(R(b,:))==L(b));
if(numel(Bb)>hb)
    Bb = Bb(1:hb);
end


%% 

%if(~iIsSingleton)
%    Ba = Ba(L(Ba)==L(i));
%end
%if(~jIsSingleton)
%    Bb = Bb(L(Bb)==L(j));
%end


Pa = min(P([a,Ba],:),[],1);
Pb = min(P([b,Bb],:),[],1);
if(isempty(Ba))
    assert(all(Pa==P(a,:)))
end
if(isempty(Bb))
    assert(all(Pb==P(b,:)))
end

% compute jaccard index
if(all(Pa==0))
    jacc = 0;
elseif(all(Pb==0))
    jacc = 0;
else
    s = sum(Pa.*Pb);
    jacc = s/(Pa*Pa'+ Pb*Pb' - s);
end

if(numel(Ba)+numel(Bb)==0)
    % assert(abs(1-jacc - Jdist(a,b))<1e-5);
end

ok = jacc> 0.01;

end