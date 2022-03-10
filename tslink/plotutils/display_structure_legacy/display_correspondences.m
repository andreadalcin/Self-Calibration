function [ ] = display_correspondences( y,G, img2, img1 )
%%
if(nargin < 4)
    imshow([img2])
    num_clust = max(G);
    cmap = brewermap(num_clust,'Set2');
    hold all;
    % map residual in [sa,sb]
    smb = 'ods';
    for i =1:max(G)
        id = G==i;
        
        u = y(1:2,id);
        v = y(4:5,id);
        for j =1: sum(id)
            line([u(1,j),v(1,j)],[u(2,j),v(2,j)],'Color', cmap(i,:),'Linewidth',0.5);
        end
        %scatter(u(1,:),u(2,:),10, cmap(i,:),'filled','MarkerEdgeColor',cmap(i,:));
        scatter(v(1,:),v(2,:),60 ,cmap(i,:),smb(mod(i,3)+1),'filled','MarkerEdgeColor','k');
        
    end
elseif(nargin == 4 )
    subplot(1,2,1);
    imshow([img1])
    num_clust = max(G);
    cmap = brewermap(num_clust,'Set2');
    hold all;
    % map residual in [sa,sb]
    smb = 'ods';
    for i =1:max(G)
        id = G==i;
        
        v = y(1:2,id);
        u = y(4:5,id);
        for j =1: sum(id)
            line([u(1,j),v(1,j)],[u(2,j),v(2,j)],'Color', cmap(i,:));
        end
        %scatter(u(1,:),u(2,:),10, cmap(i,:),'filled','MarkerEdgeColor',cmap(i,:));
        scatter(v(1,:),v(2,:),60 ,cmap(i,:),smb(mod(i,3)+1),'filled','MarkerEdgeColor','k');
        
    end
    subplot(1,2,2)
    imshow([img2])
    num_clust = max(G);
    cmap = brewermap(num_clust,'Set2');
    hold all;
    % map residual in [sa,sb]
    smb = 'ods';
    for i =1:max(G)
        id = G==i;
        
        u = y(1:2,id);
        v = y(4:5,id);
        for j =1: sum(id)
            line([u(1,j),v(1,j)],[u(2,j),v(2,j)],'Color', cmap(i,:));
        end
        %scatter(u(1,:),u(2,:),10, cmap(i,:),'filled','MarkerEdgeColor',cmap(i,:));
        scatter(v(1,:),v(2,:),60 ,cmap(i,:),smb(mod(i,3)+1),'filled','MarkerEdgeColor','k');
        
    end
end


end

