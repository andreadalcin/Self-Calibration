function [L,Z] = tlnkGric(P,X, model, gricParam, img,Y)
%TLNK  implementation with GRIC model selection
%%
if(nargin <5)
    img = [];
end

%% debug options
o_debugPlot = false;
o_verbose = false;
if(o_debugPlot)
    cmap = brewermap(2*size(P,1),'Set2');
    % get the bounding box of the data for synthetic 2D dataset
    [axlim,aylim] = get2Daxlim(X,0.01);
end
%% set merging GRIC criterion
if(strcmp(model,'line'))
    isMergeableGric = @isMergeableGricLine;
    cardmss =2;
    thCard = 3;
elseif(strcmp(model,'circle'))
    isMergeableGric = @isMergeableGricCircle;
    cardmss = 3;
    thCard = 4;
elseif(strcmp(model,'lc'))
    isMergeableGric = @isMergeableGricLC;
    cardmss = 3;
    thCard = 4;
elseif(strcmp(model,'lcp'))
    cardmss = 3;
    isMergeableGric = @isMergeableGricLCP;
    thCard = 4;
elseif(strcmp(model,'homography'))
    isMergeableGric = @isMergeableGricHomography;
    cardmss = 4;
    thCard = 4;
elseif(strcmp(model,'affine_fundamental'))
    cardmss = 4;
    isMergeableGric = @isMergeableGricAffineFundamental;
    thCard = 20;
elseif(strcmp(model,'fundamental'))
    cardmss = 8;
    isMergeableGric = @isMergeableGricFundamental;
    thCard = 9;
elseif(strcmp(model,'haf'))
    cardmss = 8;
    isMergeableGric = @isMergeableGricHAF;
    thCard = 8;
elseif(strcmp(model,'planecylinder'))
    cardmss = 3;
    isMergeableGric = @isMergeableGricPlaneCylinder;
    thCard = 5;
elseif(strcmp(model,'plane'))
    cardmss = 3;
    isMergeableGric = @isMergeableGricPlane;
    thCard = 4;
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
%% j-linkage based clustering
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
        [ok, msScore, msOutput ] = isMergeableGric(X, L, i, j, lambda1 , lambda2, sigma);
        isgric = true;
        isslnk = false;
    end
    
    %% //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    if(o_debugPlot && isgric && false)
        figure(105);
        
        clf;
        
        hold all;
        if(~isempty(img))
            displaySLMerge(Y, L,i,j, n,cmap,img);
        else
            displaySLMerge(X, L,i,j, n,cmap);
        end
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
        
        
        % update preference matrix
        P(i,:) = min(P(i,:),P(j,:));
        P(j,:) = nan;
        
        % update cluster distance matrix using the J-linkage rule
        D =  updateJlnkDist(D,P,i,j);
        % update preferences
       
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

if(o_debugPlot)
    Lpruned = prune_small_clust(L,cardmss);
    figure(1052);
    subplot(1,2,1)
    displayImageClusters(Y,Lpruned,img)
    %displayClusters(X,Lpruned);
    subplot(1,2,2);
    imagesc(P);
    pause(3);
end
end




function [dmin, row, col] = findMinDist(D,n)
% find the minimum of a square distance matrix of size n
% the diagonal of the matrix should be set to Inf
% return the minimum value, its row and col subscript indices
[dmin,index] = min(D(:));
[row,col] = ind2sub([n,n],index);
end

function D =  updateJlnkDist(D,P,i,j)
% update the cluster distance matrix
% according to the j linkage rule
for k = 1:size(D,2)
    if(all(~isnan(P(k,:))))
        D(i,k) = tanimoto_fast(P(i,:),P(k,:));
    else
         D(i,k) = Inf;
    end
end
D(:,i) = D(i,:);
D(i,i) = Inf;
D(j,:) = Inf;
D(:,j) = Inf;

end