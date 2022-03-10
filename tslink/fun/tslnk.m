function [L,Z] = tslnk(P,X)
%TSLNK  naive implementation

o_debugPlot = false;
cmap = brewermap(2*size(P,1),'Set2');
%% build tanimoto distances
n = size(P,1);
s0 = P*P'; % inner product
n0 = diag(s0); %norm
m0 = repmat(n0,1,n);
d0 = m0 + m0' - s0;
D =  1 - s0./d0;

%%
% checks on the distance matrix
if(size(D,1)~=size(D,2))
    % put the distance matrix in a squareform nxn
    D= squareform(D);
end
D = D + diag(inf.*ones(n,1)); % distance cluster

% preallocations
Z = zeros(n-1,3); % output dendrogram matrix.
L = 1:n; % clustering labelling

%% main loop
s = 0;            % iteration count
[dmin, i, j] = findMinDist(D,n);
while(dmin<1)
    % check if clusters of i and j can be merged
    if(sum(L==L(i))==1 && sum(L==L(j))==1)
        % if clusters are singleton, merge
        ok = true;
        %debugPlot = false;
    else
        [ok] = isMergeableSL( P, L, i, j);
        %debugPlot = true;
    end
    %%%%%%
    if(o_debugPlot)
        %%
        figure(105);
        clf
        displaySLMerge(X, L,i,j, n,cmap);
        pause(0.1);
      
    end
    %%%%%%
    
    if(ok)
        s = s+1;
        % update dendrogram
        Z(s,:) = [L(i) L(j) dmin];
        % update clustering
        L(L==L(i)) = n+s;
        L(L==L(j)) = n+s;
        
        % update cluster distance matrix
        D =  updateSlnkDist(D,i,j);
    else
        %fprintf('%i no merge\n',s);
        D(i,j) = Inf;
        D(j,i) = Inf;
    end
    % find minimum
    [dmin, i, j] = findMinDist(D,n);
    
    
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);
Z = Z(1:s,:); % added
L = grp2idx(L);
end




function [dmin, row, col] = findMinDist(D,n)
% find the minimum of a square distance matrix of size n
% the diagonal of the matrix should be set to Inf
% return the minimum value, its row and col subscript indices
[dmin,index] = min(D(:));
[row,col] = ind2sub([n,n],index);
end

function D =  updateSlnkDist(D,i,j)
% update the cluster distance matrix
% according to the single linkage rule
D(i,:) = min(D(i,:),D(j,:));
D(:,i) = D(i,:);
D(i,i) = Inf;
D(j,:) = Inf;
D(:,j) = Inf;

end