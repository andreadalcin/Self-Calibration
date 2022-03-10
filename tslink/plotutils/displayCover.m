function [] = displayCover(X,U,opts)
%DISPLAYCOVER display a covering in the plane
assert(~isempty(X),'Cannot display cluster of empty data');
assert(size(U,1)==size(X,2),'Dimension mismathc between number of data points and clustering vector');

%%  plotting options
if(nargin<3)
    % set default plotting options
    defaultOpts = defaultOptsDisplay();
    [mrkrSize, syms, scheme, colorOutlier] = parseOptsClustDisplay(defaultOpts);
else
    [mrkrSize, syms, scheme, colorOutlier] = parseOptsClustDisplay(opts);
end
%% chek if a figure already exists
g = groot;
assert(~isempty(g.Children),'Cannot display segmentation... have you create a figure?');

%% start the plotting
hold all;
numSeg = size(U,2); % number of segments
cmap = brewermap(numSeg,scheme); % add color for outlier

% plot outlier
outlier = sum(U,2)==0;
scatter(X(1,outlier),X(2,outlier),mrkrSize,colorOutlier,'x');

% plot inliers
for j  = 1:numSeg
    curSym = syms(mod(j-1,numel(syms))+1);
    inliers = U(:,j)>0;
    scatter(X(1,inliers),X(2,inliers),mrkrSize,cmap(j,:),curSym);
end

axis off;
axis equal;
legend off;
hold off;

end

