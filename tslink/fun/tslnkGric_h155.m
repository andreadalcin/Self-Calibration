function [L,Z] = tslnkGric_h155(P,X, model, gricParam)
%TSLNK  naive implementation


%% debug options
o_debugPlot = false;
o_verbose = false;
if(o_debugPlot)
    cmap = brewermap(2*size(P,1),'Set2');
    % get the bounding box of the data for synthetic 2D dataset
    [axlim,aylim] = get2Daxlim(X,0.01);
end
%% set is mergeable model

if(strcmp(model,'flat2'))
    cardmss = 3;
    dimFlat = 2;
    isMergeableGric = @isMergeableGricFlat;
    thCard = 3;
elseif(strcmp(model,'flat3'))
    cardmss = 4;
    dimFlat = 3;
    isMergeableGric = @isMergeableGricFlat;
    thCard = 4;
elseif(strcmp(model,'mixedFlats'))
    cardmss = 4;
    dimFlat = [2,3,4];
    isMergeableGric = @isMergeableGricFlatsMixed;
    thCard = max(dimFlat)+1;
end
% set gric parameters
lambda1 = gricParam.lambda1;
lambda2 = gricParam.lambda2;
sigma = gricParam.sigma;

%% build tanimoto distances between data
n = size(P,1);
s0 = P*P'; % inner product
n0 = diag(s0); %norm
m0 = repmat(n0,1,n);
d0 = m0 + m0' - s0;
D =  1 - s0./d0;
%% single linkage based clustering
% checks on the distance matrix
if(size(D,1)~=size(D,2))
    % put the distance matrix in a squareform nxn
    D= squareform(D);
end
D = D + diag(inf.*ones(n,1)); % D distance between clusters
D0 = D; % original distance matrix

% preallocations
Z = zeros(n-1,3); % output dendrogram matrix.
L = 1:n; % clusters labels
L = L(:);

%% main loop
s = 0;            % iteration count
[dmin, i, j] = findMinDist(D,n);
while(dmin<Inf)
    
    cardi = sum(L==L(i));
    cardj = sum(L==L(j));
    % check if clusters of i and j can be merged
    if(cardi==1 && cardj==1)
        % if clusters are singleton, merge
        ok = true;
        isgric = false;
        isslnk = false;
    elseif(cardi< thCard || cardj < thCard)
        % if clusters aren't big enough to fit a model do Single Linkage
        % test
        ok = isMergeableSL( P, L, i, j);
        isgric = false;
        isslnk = true;
    else
        % perform gric test
        [ok, msScore, msOutput ] = isMergeableGric(X, L, i, j, lambda1 , lambda2, sigma, dimFlat);
        isgric = true;
        isslnk = false;
    end
    
    %% //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    if(o_debugPlot && isgric)
        figure(105);
        
        clf;
        
        hold all;
            displaySLMerge(X, L,i,j, n,cmap);
        if(isgric)
            if(ok)
                disp(['merge accepted with ', msScore.model])
                title(['merge accepted with ', msScore.model])
               
            else
                disp(['merge rejected with ', msScore.model])
                title(['merge rejected with ', msScore.model])
            end
            drawnow;
        end
        if(isslnk)
            if(ok)
                disp('SL merge accepted');
            else
                disp('SL merge rejected');
            end
        end
       

        
        fprintf('Waiting input...\n')
        %pause;
        pause(0.01);
        
    end
    %% //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    
    if(ok)
        if(o_verbose)
            fprintf('iteration: %i\n',s);
        end
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
L = grp2idx(L);

% if(o_debugPlot)
%     Lpruned = prune_small_clust(L,cardmss);
%     figure(1052);
%     subplot(1,2,1)
% %    displayImageClusters(Y,Lpruned,img)
%     %displayClusters(X,Lpruned);
%     subplot(1,2,2);
%     imagesc(P);
%     pause(3);
% end
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