function [ ] = display_clust( X,G ,switch_out,smb,col)
%%
if(nargin<3)
    switch_out = 1;
    smb = 'o';
    col = [];
end
d = size(X,1);
if(nargin==1)
    if(d==2)
        scatter(X(1,:),X(2,:),'.')
    elseif(d==3)
        scatter3(X(1,:),X(2,:),X(3,:),'.')
    end
else
    % if a segmentation is provided
    num_clust = max(G);
    cmap = brewermap(num_clust,'Set2');
    hold all;
    
    for i =1:max(G)
        id = G==i;
        if(d==2) %2d plotting
            x = X(1,id);
            y = X(2,id);
            if(isempty(col))
               % scatter(x,y,80,cmap(i,:),'filled','Marker',smb,'MarkerEdgeColor',cmap(i,:));
               scatter(x,y,80,cmap(i,:),'filled','Marker',smb,'MarkerEdgeColor','k');
            else
                scatter(x,y,50,col,'Marker',smb,'MarkerEdgeColor',col);
            end
            if(switch_out==1)
                scatter(X(1,G==0),X(2,G==0),50,[0.3,0.3,0.3],'filled','Marker',smb,'MarkerEdgeColor',[0.2,0.2,0.2]);
            end
            
            % if you want label numbers
            %name = num2str(i.*ones(numel(x),1));
            %dx = 0.01; dy = 0.01; % displacement so the text does not overlay the data points
            %text(x+dx, y+dy, name, 'Color',cmap(i,:));
        elseif(d==3) %3d plotting
            x = X(1,id);
            y = X(2,id);
            z = X(3,id);
            if(isempty(col))
                scatter3(x,y,z,20,cmap(i,:),'filled','MarkerEdgeColor',cmap(i,:));
            else
                scatter3(x,y,z,20,col,'filled','MarkerEdgeColor',col);
            end
           % name = num2str(i.*ones(numel(x),1));
           % dx = 0.01; dy = 0.01; % displacement so the text does not overlay the data points
            %text(x+dx, y+dy,z, name, 'Color',cmap(i,:));
        elseif(d==6)
            x = X(1,id);
            y = X(2,id);
            scatter(x,y,10,cmap(i,:),'filled','Marker',smb,'MarkerEdgeColor',cmap(i,:));
        end
        
    end
end

