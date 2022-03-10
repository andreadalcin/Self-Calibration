function x = hash(i,n,p)

% famiglia paremetrica di funzioni di hash, con parametro i
% restituisce una approssimazione di una permutazione random di [0... n-1]
% in realta' mappa [0... n-1] in [0... p-1]
% p e' il prossimo intero maggiore di n; viene passato come parametro solo
% per efficienza

% %  esempio con due funzioni di has
%  h1 =@(x) mod(x+1,5);
% % 
%  h2 =@(x) mod(3*x+1,5);
% 
% % 
%  if i==1
%      x = h1([0:n-1]);
%  else
%      x = h2([0:n-1]);
% % 
%  end


 h1 =@(x) mod(2*x+3,p);
 
 h2 =@(x) mod(5*x+7,p);

  x =  mod(h1([0:n-1]) +  i * h2([0:n-1]), p);
 
