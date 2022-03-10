function [ C, card ] = keep_k_clust(C, kappa )
%PRUNE_SMALL_CLUST prune small clusters: set to zero segments with
% theta points or less
 lbl = sort(unique(C));
 card = numel(lbl);
 
 for i = 1:numel(lbl)
    id = (C==lbl(i));
    card(i) = sum(id);        
 end
 
[card ,o]= sort(card,'descend');
 
for i = kappa+1:numel(lbl)
  C(C==o(i)) = 0; 
end

C(C~=0)=grp2idx(C(C~=0));
 
end

