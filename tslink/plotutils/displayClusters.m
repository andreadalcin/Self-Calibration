function [] = displayClusters(X,F,opts)
%DISPLAYCLUSTERS display a clustering in the plane
assert(size(F,1)==size(X,2),'Dimension mismatch between number of data points and clustering vector');
if(isempty(X))
    return;
end

assert(~isempty(X),'Cannot display cluster of empty data');


%%  plotting options
if(nargin<3)
    % set default plotting options
    defaultOpts = defaultOptsClustDisplay();
    [mrkrSize, syms, scheme, colorOutlier,~,mrkrFaceAlpha] = parseOptsClustDisplay(defaultOpts);
else
    [mrkrSize, syms, scheme, colorOutlier,~,mrkrFaceAlpha] = parseOptsClustDisplay(opts);
end
%syms = 'o';
scheme = 'Set2';
mrkOutlier = 'o';
%% chek if a figure already exists
g = groot;
assert(~isempty(g.Children),'Cannot display segmentation... have you created a figure?');

%% perform the plotting
hold all;
numSeg = max(F); % number of segments
cmap = brewermap(numSeg,scheme); % add color for outlier

% plot outlier
d = size(X,1);
if(d==2)
    scatter(X(1,F==0),X(2,F==0),mrkrSize,colorOutlier,mrkOutlier,'filled');
elseif(d>=3)
    scatter3(X(1,F==0),X(2,F==0),X(3,F==0),mrkrSize,colorOutlier,mrkOutlier);
end
% plot inliers
for j  = 1:numSeg
    curSym = syms(mod(j-1,numel(syms))+1);
    if(mod(j-1,numel(syms)==0))
        syms = syms(randperm(numel(syms)));
    end
    if(d==2)
        scatter(X(1,F==j),X(2,F==j),mrkrSize,cmap(j,:),curSym,'filled','MarkerEdgeColor','k','LineWidth',0.5,'MarkerFaceAlpha',mrkrFaceAlpha);
    elseif(d>2)
        scatter3(X(1,F==j),X(2,F==j),X(3,F==j),mrkrSize,cmap(j,:),curSym);
    end
end

%axis off;
%axis equal;
%legend off;
%hold off;

end


