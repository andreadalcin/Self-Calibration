function [ ] = display_lines_flg( X, F, M,epsi,flg )

hold all;
col_map = brewermap(max(F),'Set2');
for i =1:max(F)
    if(flg(i))
        if(sum(F==i)>0)
            display_band(X(:,F==i), M(:,i),epsi,col_map(i,:));
            scatter(X(1,F==i),X(2,F==i),5,col_map(i,:),'filled','MarkerEdgeColor',col_map(i,:))
        end
    end
    hold on
end
%display_clust(X,F);



end

