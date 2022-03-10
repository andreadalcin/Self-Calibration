function [X,G,x] = load_data_h155( sequencename, out )
%LOAD_DATA_H155 Summary of this function goes here
%  out=1 load outliers
%  out=0 load inliers only

if (nargin<2)
    % load inlier data
    out = 0;
end

    load(strcat(sequencename,'_truth.mat'));
if(~exist('frames','var'))
    frames = size(x,3);
end
if(~exist('points','var'))
    points = size(x,2);
end
X = [];
for f= 1:frames
    X = vertcat(X, x(1:2,1:points,f));
end
%%

%load('pathH155.mat')

%%

if(out==1)
    [xoutliers,youtliers] = generate_outliers_randomwalk2(sequencename);
    %load(strcat(sequencename,'/outliers_data.mat'));
    Xout=[];
    for f= 1:frames
        Xout = vertcat(Xout, xoutliers(1:2,1:points,f));
    end
    X=[X,Xout(:,[1:ceil(points*0.2)])];
end

%save(strcat(sequencename,'/outliers_data.mat'));
% prepare ground-truth
G=zeros(size(X,2),1);
G(1:points)=s';


% reorder data
[~, order]=sort(G);
G=G(order);
X=X(:,order);



end

