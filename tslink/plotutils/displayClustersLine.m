function [] = displayClustersLine(X,F,opts, epsi)
%DISPLAYCLUSTERS display a clustering in the plane
assert(size(F,1)==size(X,2),'Dimension mismatch between number of data points and clustering vector');

if(isempty(X))
    return;
end
if(nargin< 4)
    % automatically set the bandwidht of the inlier band to the maximum
    % residual of the model
    autoEpsi = true;
else
    autoEpsi = false;
end
assert(~isempty(X),'Cannot display cluster of empty data');
F = prune_small_clust(F,2);
M = recover_line(X,F);

%%  plotting options
if(isempty(opts))
    % set default plotting options
    defaultOpts = defaultOptsClustDisplay();
    [mrkrSize, syms, scheme, colorOutlier,~,mrkrFaceAlpha] = parseOptsClustDisplay(defaultOpts);
else
    [mrkrSize, syms, scheme, colorOutlier,~,mrkrFaceAlpha] = parseOptsClustDisplay(opts);
end

%% chek if a figure already exists
g = groot;
assert(~isempty(g.Children),'Cannot display segmentation... have you create a figure?');

%% perform the plotting
hold all;
numSeg = max(F); % number of segments
cmap = brewermap(numSeg,scheme); % add color for outlier

% plot outlier
d = size(X,1);

scatter(X(1,F==0),X(2,F==0),mrkrSize,colorOutlier,'x');
% plot bands
for j  = 1:numSeg
    curSym = syms(mod(j-1,numel(syms))+1);
    if(mod(j-1,numel(syms)==0))
        syms = syms(randperm(numel(syms)));
    end
    if(autoEpsi)
        resi = res_line(X(:,F==j), M(:,j));
        bandWidth = x84_(resi,1);
        bandWidth = max(resi);
    else
        bandWidth = epsi;
    end
    display_band(X(:,F==j), M(:,j),bandWidth,cmap(j,:));
end

% plot inliers
for j  = 1:numSeg
    curSym = syms(mod(j-1,numel(syms))+1);
    if(mod(j-1,numel(syms)==0))
        syms = syms(randperm(numel(syms)));
    end
    scatter(X(1,F==j),X(2,F==j),mrkrSize,cmap(j,:),curSym,'filled','MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',mrkrFaceAlpha);
end

axis off;
axis equal;
legend off;
hold off;

end


