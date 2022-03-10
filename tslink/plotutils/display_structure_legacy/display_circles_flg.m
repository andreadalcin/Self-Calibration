function [ ] = display_circles_flg( X,F,M,epsi,flg )

scatter(X(1,:),X(2,:),nan,'.','w')
hold on;
axis equal;
ax = gca;
width_x = ax.XLim(2)-ax.XLim(1);
width_y = ax.YLim(2)-ax.YLim(1);
xlim([ax.XLim(1) - width_x,ax.XLim(2)+width_x]);
xlim([ax.YLim(1) - width_y,ax.YLim(2)+width_y]);
axis manual
col_map = brewermap(max(F),'Set2');
for i =1:max(F)
    if(flg(i))
    display_anulus(X(:,F==i),M(:,i),epsi,col_map(i,:));
     scatter(X(1,F==i),X(2,F==i),5,col_map(i,:),'filled','MarkerEdgeColor',col_map(i,:));
    hold on
    end
end





end

