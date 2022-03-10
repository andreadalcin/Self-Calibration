function [ ] = display_homography_residual( X,G,y, img1, img2, res_gt )
%%

imshow([img1,img2])
num_clust = max(G);
cmap = brewermap(num_clust,'Set2');

hold all;
[~,res] = recover_homography( X,G );
if(nargin<6)
    res_gt = res;
end
% map residual in [sa,sb]
sa = 50;
sb = 300;
ra = min(res_gt);
rb = max(res_gt);
ms = (sb-sa)/(rb-ra);
s = ms.*res + sa - ms*ra;

ta = 0.9;
tb = 1;
ta = min(res);
tb = max(res);
mt = (tb-ta)/(rb-ra);
t = mt.*res + ta - mt*ta;

t(isnan(t))= 1;

for i =1:max(G)
    id = G==i;
    
    u = y(1,id);
    v = y(2,id);
    scatter(u,v,s(id),cmap(i,:),'filled','Marker','s','MarkerEdgeColor','k');
    u = y(4,id)+size(img1,2);
    v = y(5,id);
    scatter(u,v,s(id),cmap(i,:),'filled','Marker','s','MarkerEdgeColor','k');
 
end


end

