function [] = displaySynthHomo(X, E, G)
%DISPLAYSYNTHHOMO
mrkrSize = 40;
if(nargin ==1)
    displayGrid(X(1:2,:), E); 
    displayGrid(X(4:5,:), E); 
else
    cmap = brewermap(max(G),'Accent');
    hold all;
    for j = 1:max(G)
        clr = cmap(j,:);
        % plot grid
        % select only the edge that joining vertex belonging to the cluster
        Eclust = E(G(E(:,1))==j & G(E(:,2))==j,:);
        displayGrid(X(1:2,:), Eclust, clr); 
        displayGrid(X(4:5,:), Eclust, clr); 
        % plot points
        scatter(X(1,G==j),X(2,G==j),mrkrSize,'o','MarkerEdgeColor',clr,'MarkerFaceColor',clr);
        scatter(X(4,G==j),X(5,G==j),mrkrSize,'o','MarkerEdgeColor',clr,'MarkerFaceColor',clr);
    end
    scatter(X(1,G==0),X(2,G==0),'kx');
    scatter(X(4,G==0),X(5,G==0),'kx');
end
axis equal;
end

