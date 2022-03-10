function labelHp = getMssLabels(S,F)
%GETMSSLABELS: given a MSS-sampling of the hypotheses and a clustering vector F
% this function gives a label to each MSS according to F.
% More specifically, each mss can be
% (a) pure: all the sample are taken from the same inlier cluster and get the
%   label of that cluster.
% (b)impure:
%    (b1) there is at least an outlier, the mss is labelled as 0
%    (b2) all the sample are inliers but from different clusters, labelled as -1


m = size(S.S,2);
labelHp = -1.*ones(m,1);
for j = 1:m
    if(all(F(S.S(2:end,j)) == F(S.S(1,j))))
        labelHp(j) = F(S.S(1,j));
    elseif(any(F(S.S(:,j))==0))
        labelHp(j) = 0;
    end
    
end
end

