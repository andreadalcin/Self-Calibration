function [ok,Ba,Bb,a,b] = isMergeable(D, P, L, i, j, radius)
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
%   - P: the preference matrix
%   - L: vector of clusters labels
%   - i and j: index of points belonging to the cluster to be merged
%   - radius: the radius of the jaccard ball
%   - idx: auxiliary variable of points idx
%  Output:
%   - ok: tells if the cluster can be merged 


isInCi = L==L(i);
isInCj = L==L(j);
% find closest pair (a,b)
iIsSingleton = sum(isInCi)==1;
jIsSingleton = sum(isInCj)==1;
if(~iIsSingleton || ~jIsSingleton)
   [a,b] = findClosestInClusters(D, isInCi,isInCj);
else
    a = i;
    b = j;
end
assert(L(a)==L(i));
assert(L(b)==L(j))

% find points in ball centered in the closest pair
Ba = getIdsBall(a ,radius, D);
Bb = getIdsBall(b ,radius, D);
%k = 5;
%Ba = getIdsNN(a ,k, Jdist);
%Bb = getIdsNN(b ,k, Jdist);
Ba = Ba(L(Ba)==L(i));
Bb = Bb(L(Bb)==L(j));

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

ok = jacc> 0.0;

end