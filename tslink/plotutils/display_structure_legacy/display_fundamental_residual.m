function [ ] = display_fundamental_residual( X,G,y, img1, img2 )
%%
smb = 'ods';
d = size(X,1);
imshow([img1,img2])
num_clust = max(G);
cmap = brewermap(num_clust,'Set3');
hold all;
[~,res] = recover_fundamental( X,G );
% map residual in [sa,sb]
sa = 50;
sb = 300;
ra = min(res);
rb = max(res);
ms = (sb-sa)/(rb-ra);
s = ms.*res + sa - ms*ra;

for i =1:max(G)
    id = G==i;
    
    u1 = y(1,id);
    v1 = y(2,id);
    scatter(u1,v1,s(id),cmap(i,:),smb(mod(i,3)+1),'filled','MarkerEdgeColor','k');
     text(u1,v1,num2str(i),'Color','White','fontsize',5);
    u2 = y(4,id)+size(img1,2);
    v2 = y(5,id);
    scatter(u2,v2,s(id),cmap(i,:),smb(mod(i,3)+1),'filled','MarkerEdgeColor','k');
    text(u2,v2,num2str(i),'Color','White','fontsize',5);
    %for j=1:numel(u1)
        %line([u1(j),u2(j)],[v1(j),v2(j)])
    %end
end


end

