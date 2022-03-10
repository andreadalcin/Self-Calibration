function [Fout] = adjustIntersectionNeigh(X,F,R,epsi, kNeigh)
%ADJUSTINTERSECTION assigns to each inlier in an intersection the most common label in that intersection.
% every time an inlier is in the intersection of two or more models
% it is assigned to the model of the majority of its neighborhood.


% Costruisce un kdtree...

debugPlot = true;

n = size(X,2);
kappa = size(R,2);
idxModel = 1:kappa;
if(nargin < 5)
    kNeigh = 3; % number of neighbouring point searched
end
Fout = F;

tree = KDTreeSearcher(X');

for i = 1:n
    if(F(i)==0)
        % skip outlier
        continue;
    end
    isSharedModel = R(i,:)<epsi; % models that explain the point
    idxSharedModel = idxModel(isSharedModel); % idx of the model the point belongs to
    if(sum(isSharedModel)>1)
        % Find the labels of neighbouring points
        idxNeigh = knnsearch(tree,X(:,i)','k',kNeigh);
        neighLabels = F(idxNeigh);
        % Consider only the labels of the models that explain the point
        isSharedLabels =  ismember(neighLabels,idxSharedModel);
        labels = neighLabels(isSharedLabels);
        Fout(i) = mode(labels);
        if(debugPlot)
            if(F(i)~=Fout(i))
                figure;
                gscatter(X(1,:),X(2,:),F);
                hold all;
                plot(X(1,i),X(2,i),'k+','MarkerSize',50);
                scatter(X(1,idxNeigh),X(2,idxNeigh));
                legend off;
            end
        end
    end
    
    
end


end

