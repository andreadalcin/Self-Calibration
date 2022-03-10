function [L,Z] = tslnkGric(P,X, model, gricParam, img,Y)
%TSLNK  naive implementation
%%
if(nargin <5)
    img = [];
end

%% debug options
o_debugPlot = true;
o_verbose = false;
o_save_gif = true;
if(o_debugPlot)
    cmap = brewermap(2*size(P,1),'Set2');
    % get the bounding box of the data for synthetic 2D dataset
    [axlim,aylim] = get2Daxlim(X,0.01);
    filename = 'clustering.gif';
    count = 1;
    
end
%% set is mergeable model
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
    thCard = 3;
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
    thCard = 3;
elseif(strcmp(model,'sphere'))
    cardmss = 4;
    isMergeableGric = @isMergeableGricSphere;
    thCard = 4;
elseif(strcmp(model,'planecylindersphere'))
    cardmss = 4;
    isMergeableGric = @isMergeableGricPlaneCylinderSphere;
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
        [ok, msScore, msOutput ] = isMergeableGric(X, L, i, j, lambda1 , lambda2, sigma);
        isgric = true;
        isslnk = false;
    end
    
    %% //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
    if(o_debugPlot && isgric)
        %if(o_debugPlot && isgric && s>50)
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        clf;
        colorAccept ='g';
        colorReject = 'r';
        hold all;
        
        if(isgric)
            if(ok)
                
                disp(['merge by ', msScore.model])
                title(['merge clusters by ', msScore.model],'fontsize',30,'fontname','Source Serif Pro')
                if(strcmp(msScore.model,'line'))
                    display_band( X,msOutput.mij(:) , sigma,colorAccept,0.7)
                elseif(strcmp(msScore.model,'circle'))
                    display_anulus( X, msOutput.mij(:), sigma,colorAccept,0.7)
                end
            else
                disp(['merge rejected by ', msScore.model])
                title(['\color{red}reject merge'],'fontsize',30,'fontname','Source Serif Pro')       
                
            end
            
        end
        bb = getBB(X,0.1);
        set(h,'color','white');
        
        
        
        if(~isempty(img))
            displaySLMerge(Y, L,i,j, n,cmap,img);
        else
            displaySLMerge(X, L,i,j, n,cmap);
        end
        xlim([bb.xmin,bb.xmax]);
        ylim([bb.ymin,bb.ymax]);
        drawnow;
        if(isslnk)
            if(ok)
                disp('SL merge accepted');
            else
                disp('SL merge rejected');
                pause(1)
            end
        end
        %         if(~isempty(img))
        %
        %             subplot(2,2,1)
        %             imshow(img); hold all;
        %             if(~ok)
        %                 scatter(Y(1,msOutput.isInCi),Y(2,msOutput.isInCi),200,'rs','filled');
        %                 scatter(Y(1,msOutput.isInCj),Y(2,msOutput.isInCj),200,'mo','filled');
        %             else
        %                 scatter(Y(1,msOutput.isInCi),Y(2,msOutput.isInCi),200,'gs','filled');
        %                 scatter(Y(1,msOutput.isInCj),Y(2,msOutput.isInCj),200,'co','filled');
        %             end
        %             subplot(2,2,3);
        %             displayImageClusters(Y,L,img);
        %             subplot(2,2,[2,4])
        %             imagesc(L);
        %
        %         else
        %             scatter(msOutput.Xi(1,:),msOutput.Xi(2,:),200,'rs');
        %             scatter(msOutput.Xj(1,:),msOutput.Xj(2,:),200,'mo');
        %             if(strcmp(msScore.model,'line'))
        %                 if(ok)
        %                     display_band( msOutput.Xij, msOutput.mij, max(msOutput.rij),'g' );
        %                 else
        %                     display_band( msOutput.Xi, msOutput.mi, max(msOutput.ri),'r' );
        %                     display_band( msOutput.Xj, msOutput.mj, max(msOutput.rj),'m' );
        %                 end
        %             elseif(strcmp(msScore.model,'circle'))
        %                 if(ok)
        %                     displayAnularBand( msOutput.Xij, msOutput.mij, max(msOutput.rij),'g' );
        %                 else
        %                     displayAnularBand( msOutput.Xi, msOutput.mi, max(msOutput.ri),'r' );
        %                     displayAnularBand( msOutput.Xj, msOutput.mj, max(msOutput.rj),'m' );
        %                 end
        %             end
        %             scatter(X(1,:),X(2,:),'k.');
        %             xlim(axlim);
        %             ylim(aylim);
        %         end
        
        
        
        
        fprintf('Waiting input...\n')
        %pause;
        set(gcf,'color','w');
        %pause(0.5);
        
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
    if(o_save_gif && isgric)
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if count == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            count = count+1;
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);
L = grp2idx(L);

if(o_debugPlot)
    %     Lpruned = prune_small_clust(L,cardmss);
    %     figure(1052);
    %     subplot(1,2,1)
    %     displayImageClusters(Y,Lpruned,img)
    %     %displayClusters(X,Lpruned);
    %     subplot(1,2,2);
    %     imagesc(P);
    %     pause(3);
    
    Lpruned = prune_small_clust(L,cardmss);
   h = figure('units','normalized','outerposition',[0 0 1 1]);
    
    clf;
    optsClustSpy = defaultOptsClustDisplay();
    optsClustSpy.syms = 'o';
    optsClustSpy.colorOutlier = [0.2,0.2,0.2];
    optsClustSpy.mrkrSize = 100;
    displayClusters(X,Lpruned,optsClustSpy);
    axis equal;
    axis off;
    set(gcf,'color','white')
    xlim([bb.xmin,bb.xmax]);
    ylim([bb.ymin,bb.ymax]);
    frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');
end
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