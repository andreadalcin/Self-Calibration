function [C,T,P,active] = jlnk_minhash(P,k)
%JLNK perform j-linkage clustering on the preferenece matrix P
% approximating jaccard distances with minhash
% author: Luca Magri
% date march 2018
% improoved efficiency thanks to nearest neighbors pairs
% this implementation was inspired by
% - Fast and Memory Efficient Implementationof the Exact PNN
% - Modern hierarchical, agglomerative clustering algorithms
% - Fast Agglomerative Clustering Using a k-Nearest Neighbor Graph

%%
[n,m] = size(P);      % number of points
active = true(n,1); % active labels of cluster
T = nan(n-1,3);     % dendrogram
C = (1:n)';
%%
if(nargin==1)
    k = min(400,m);
end
%%
p = nextprime(m); % il prossimo primo superiore a n
Pmh = zeros(n,k); % Preference matrix versione minhash
H = nan(m,k);
for  j = 1:k
    h = hash(j,m,p); % per ogni j approssima una random permutation di [0...m]
    H(:,j) = h;
    for i = 1:n
        
        Pmh(i,j) = min(h(P(i,:)==1)); % minhash signature
        
    end
end
%% compute nearest neighbors
% nnD  - distance from a point to its nearest neighbor
% nnID - id of the nearest neighbor
% E    - nearest pairs: distance - cluster id - nn id
%D = squareform( pdist(P,'jaccard'));
D = squareform( pdist(Pmh,'hamming'));
D = D + diag(inf.*ones(n,1)); % remove self distances
[nnD,nnId]= min(D,[],2);
E = [nnD(:,1), (1:n)', nnId(:,1)];
E = sortrows(E,1); % keep the minimum distance on top

%% main loop
niter = 0;
d_lazy = inf; % keep the minimum distance with new created clusters for lazy sorting of E
converged = false;

while(~converged)
    
    niter = niter+1;
    
    
    %% find closest cluster
    
    x = E(1,2); % id of the first cluster and id of the new merged cluster
    y = E(1,3); % id of the second cluster (it will become unactive)
    d = E(1,1); % minimum distance
    
    if(d == 1)
        break
    end
    
    if(x>y) % maintain x<y
        temp = x;
        x = y;
        y = temp;
    end
    T(niter,:) = E(1,:); % happend to dendrogram
    
    %% merge pair
    
    P(x,:) = P(x,:).*P(y,:);           % update preference vector of new cluster
    P(y,:) = P(x,:);                   % optional
    
    % hashed representation
    for j=1:k
        h = H(:,j);
        
        Pmh(x,j) = min(h(P(x,:)==1));
        
    end
    Pmh(y,:) = Pmh(x,:);
    
    active(y) = false;                 % deactivate y
    
    
    C(C==x| C==y) = x;
    %% compute new distances
    
    % compute distances between xUy and other active cluster
    % and store the minimum distance in dmin
    % update also d_lazy
    
    dmin = inf;
    for i = 1:n
        if(active(i) && i~=x)
            %d_xi = 1 - sum(P(x,:).*P(i,:))/sum(P(x,:)+P(i,:)>0); % jaccard distance
            %d_xi1 = pdist2(Pmh(x,:),Pmh(i,:),'hamming');
            d_xi = sum(Pmh(x,:)~=Pmh(i,:))/k;
            D(x,i) = d_xi;
            
            if(d_xi < dmin)
                dmin  = d_xi;
                neigh = i;
                if(dmin < d_lazy)
                    d_lazy = dmin;
                end
            end
        else
            D(:,i) = inf;
        end
    end
    D(:,x) = D(x,:);
    
    E(E(:,2)==x | E(:,2)==y ,:) = [];  % remove x and y from nearest pair
    E(end+1,:) = [dmin, x, neigh];     % add new merged cluster and label it as x
    
    
    %% update nearest neighbor pairs
    
    for i = 1:size(E,1)-1
        id = E(i,2); % cluster id
        
        assert(all(isinf(D(id,~active))))
        
        if(~active(id))
            keyboard
        end
        if(E(i,3)== x || E(i,3)== y)           % se x o y erano vicini a id
            
            if(D(x,id)<= E(i,1))               %   if xUy  si e' avvicinato
                E(i,1) = D(x,id);              %   aggiorno la distanza
                E(i,3) = x;                    %   uso l'etichetta x (y rimossa)
            else                               %   altrimenti se xUy si e' allontanato
                [E(i,1),E(i,3)]= min(D(id,:)); %   cerco un nuovo vicino
            end
            
            
        else                                   % altrimenti se il vicino di id non ? x o y
            
            if(D(x,id)<= E(i,1))               % controllo se xUy e' piu' vicino
                E(i,1) = D(x,id);              % in tal caso aggiorno la distanza
                E(i,3) = x;                    % e l'etichetta
            end
            
        end
    end
    
    %% sort rows
    if(d_lazy < E(1,1))
        E = sortrows(E,1);
        d_lazy = inf;
    end
    
    converged  = (d==1);
end


%% extract labeling from the dendrogram

C = grp2idx(C);




end

