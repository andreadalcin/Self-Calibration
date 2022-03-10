function [ ] = display_lines( X, F, M,epsi )

hold all;
col_map = brewermap(max(F),'Set2');
for i =1:max(F)
    if(sum(F==i)>0)
        display_band(X(:,F==i), M(:,i),epsi,col_map(i,:));
       
    end
    hold on
end
display_clust(X,F);



end

