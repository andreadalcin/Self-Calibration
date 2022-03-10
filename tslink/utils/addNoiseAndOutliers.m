function [Xperturbed] = addNoiseAndOutliers(X, sigma,fracOutlier)
%ADDNOISEANDOUTLIERS take a clean dataset and add noise and outliers
% add noise
X = X + sigma.*randn(size(X));
% add outliers
nInliers = size(X,2);
nOutliers = ceil(fracOutlier*nInliers/(1-fracOutlier));

k = 0.1; % expand the bounding box
if(size(X,1)==2)
xmin = min(X(1,:));
xmax = max(X(1,:));
ymin = min(X(2,:));
ymax = max(X(2,:));
wx = xmax-xmin;
wy = ymax -ymin;
dx = 0.1*wx;
dy = 0.1*wy;
X = [X,[(xmax-xmin+2*dx ).*rand(1,nOutliers)+(xmin-dx); (ymax- ymin + 2*dy).*rand(1,nOutliers)+(ymin -dy)]];
G = [G; zeros(nOutliers,1)];
n = size(X,2); % number of points;
end

