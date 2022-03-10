function [F] = prune_outlier_cardgap(C,cardmss)
k = max(C);
card = zeros(k+1,1);
gap  =zeros(k,1);
for i = 1:k
    card(i) = sum(C==i);
end
card(k+1) = cardmss;
[card,o] = sort(card,'descend');
for i = 1:k
gap(i) = card(i) -card(i+1);
end
[maxGap, pos] = max(gap);
F = C;
if(pos ==k+1)
  return;
end
for i = pos+1:k
    F(F==o(i))=0; 
end
F(F~=0)=grp2idx(F(F~=0));