function [Fpruned] = keepKappaStructures(F,kappa)
%KEEPBIGGERSTRUCTURES outlier rejection: keep the kappa bigger structures
numClust = max(F);
if(numClust<kappa)
    Fpruned = F;
    return;
end
card = zeros(1,numClust);
for i = 1:numClust
    card(i) = sum(F==i);
end
Fpruned = zeros(size(F));
[~, labels] = sort(card,'descend');
for i = 1:kappa
    Fpruned(F==labels(i)) = labels(i);
end

Fpruned(Fpruned>0) = grp2idx(Fpruned(Fpruned>0));
end

