function [ miss ,  fnr, acc, gtgamma] = approxMissErrorFast( c, gt )
%FUCKFUCKME compute the me, the right way
% just skip if there are too many models!


k = max(gt); % number of ground-truth clusters
%keep only the biggest n clusters in c;
c = keep_k_clust(c,k);
[miss,gtgamma]=computeMissRatePrinci(c,gt);
fnr = nan;
acc = nan;
end

