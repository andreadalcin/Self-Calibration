function [] = display_cylinders(X,F,M)
%DISPLAY_CYLINDER 
hold all;
col_map = brewermap(max(F),'Set2');
for i =1:max(F)
    if(sum(F==i)>2)
        model = cylinder_model(M(:,i));
        plot(model)
        h = findobj(gca,'Type','Surface');
        set(h(1),'FaceAlpha',0.3,'EdgeColor',col_map(i,:),'FaceColor',col_map(i,:))
        d = res_pc_cylinder(X(:,F==i),M(:,i));
        %scatter3(X(1,F==i),X(2,F==i),X(3,F==i),50,100*d)
        %max(d)
        %pause
    end
    hold on
end
display_clust(X(1:3,:),F);






%h = findobj(gca,'Type','Line');

end

