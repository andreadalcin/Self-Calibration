function [ ] = display_circles( X,F,M,epsi )

scatter(X(1,:),X(2,:),nan,'.','w');
axis equal;
ax = gca;
width_x = ax.XLim(2)-ax.XLim(1);
width_y = ax.YLim(2)-ax.YLim(1);
xlim([ax.XLim(1) - width_x,ax.XLim(2)+width_x]);
xlim([ax.YLim(1) - width_y,ax.YLim(2)+width_y]);
axis manual
hold all;

col_map = brewermap(max(F),'Accent');
for i =1:max(F)
    display_anulus(X(:,F==i),M(:,i),epsi,col_map(i,:));
    hold on
end
display_clust(X,F);




end

