function [] = displayPrefNeigh(X,dists,radiusNeigh,radiusEpsi,cmap)
%DisplayPrefNeigh given a distance vector show 
assert(size(X,2)== numel(dists),'Dimension mismatch: for each point a distance must be specified.');
thetaResolution = 2;
theta=(0:thetaResolution:360)'*pi/180;
costheta = cos(theta);
sintheta = sin(theta);
for i = 1:size(X,2)
    value = dists(i);
    if(value<=radiusNeigh)
        col = valueToColor(value,cmap);
        displayCircPatch(X(:,i),radiusEpsi,col,costheta, sintheta);
    end
end
colormap(cmap);
colorbar;
end

function displayCircPatch(center, radius, col, costheta, sintheta)
if(nargin <4)
    thetaResolution = 2;
    theta=(0:thetaResolution:360)'*pi/180;
    costheta = cos(theta);
    sintheta = sin(theta);
end


x = center(1)+ radius.*costheta;
y = center(2)+ radius.*sintheta;

patch(x,y,col,'FaceAlpha',0.5,'EdgeAlpha',0.0) ;
axis equal
end

function col = valueToColor(value,cmap)
%minValue = 0;
%maxValue = 1;
% because we are dealing with Jaccard and Tanimoto distances
l = size(cmap,1);
step = 1/l;
id = floor(value/step);
if(id ==0)
    id = 1;
end
col = cmap(id,:);

end