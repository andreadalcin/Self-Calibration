function [] = displayImageClusters(X,F,img,opts)
%DISPLAYIMAGECLUSTERS show clusters on a single image


%% Parsing
if(nargin<4)
    % set default plotting options
    defaultOpts = defaultOptsClustDisplay();
    [mrkrSize, syms, scheme, colorOutlier, imageAlpha] = parseOptsClustDisplay(defaultOpts);
else
    [mrkrSize, syms, scheme, colorOutlier,imageAlpha] = parseOptsClustDisplay(opts);
end

%% chek if a figure already exists
g = groot;
assert(~isempty(g.Children),'Cannot display segmentation... have you create a figure?');

%% perform the plotting
hold all;
numSeg = max(F); % number of segments
cmap = brewermap(numSeg,scheme); % add color for outlier
imshow(img);
alpha(imageAlpha);
hold all;
% plot outlier
scatter(X(1,F==0),X(2,F==0),mrkrSize,colorOutlier,'x');

% plot inliers
for j  = 1:numSeg
    curSym = syms(mod(j-1,numel(syms))+1);
    scatter(X(1,F==j),X(2,F==j),mrkrSize,cmap(j,:),curSym,'filled');
end

axis off;
hold off;

end

