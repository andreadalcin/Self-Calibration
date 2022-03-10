function [ ] = display_epi_correspondences( y,G, img2,M )
%%




imshow([img2])

num_clust = max(G);
cmap = brewermap(num_clust,'Set2');
hold all;
% map residual in [sa,sb]

for i =1:max(G)
    
    F = reshape(M(:,i),[3,3]);
    epiLines = epipolarLine(F,y(1:2,:)');
    points = lineToBorderPoints(epiLines,size(img2));
    line(points(:,[1,3])',points(:,[2,4])','Color', cmap(i,:));
    
    
    id = G==i;
    
    %u = y(1:2,id);
    v = y(4:5,id);
    %scatter(u(1,:),u(2,:),10, cmap(i,:),'MarkerEdgeColor',cmap(i,:));
    scatter(v(1,:),v(2,:),20 ,cmap(i,:),'filled','MarkerEdgeColor','k');
 
end


end

