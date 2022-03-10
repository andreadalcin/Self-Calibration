function [] = displayMerge(X, L,a,b, balli, ballj, n,cmap)
%DISPLAYMERGE display the merge during the jslinkage clustering
mrkSize = 45;
mrkSizeBig = 1.3* mrkSize;
mrkSizePair = 95;
%scheme = 'Set2';
syms = 'dphv<^>';
%% 
set(gcf,'color','w');
axis equal;
axis off;
hold all;
% set drawing options for the merging pairs
%cmap = brewermap(max(L),scheme);
if(L(a)<n) % not yet clusterized
    c1 ='k';
else
    c1 = cmap(L(a),:);
end
if(L(b)<n) % not yet clusterized
    c2 ='k';
else
    c2 = cmap(L(b),:);
end
% plot balls
if(~isempty(balli))
    scatter(X(1,balli),X(2,balli),250,'o','MarkerFaceColor',c1,'LineWidth',eps,'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',0);
end
if(~isempty(ballj))
    scatter(X(1,ballj),X(2,ballj),250,'o','MarkerFaceColor',c2,'LineWidth',eps,'MarkerFaceAlpha',.7,'MarkerEdgeAlpha',0);
end
% plot not clusterized points
scatter(X(1,L<=n),X(2,L<=n),20,[0,0,0],'filled'); legend off;
% plot clusterized points
for u  = n+1:max(L)
    if(u ~=L(a) && u~=L(b))
       curSym = syms(mod(u-1,numel(syms))+1);
     scatter(X(1,L==u),X(2,L==u),mrkSize,cmap(u,:),curSym,'filled','MarkerEdgeColor',[0.2,0.2,0.2],'LineWidth',0.1,'MarkerFaceAlpha',0.8);
    end
end

% plot current clusters
scatter(X(1,L==L(a)),X(2,L==L(a)),mrkSizeBig,c1,'o','filled','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceAlpha',1);
scatter(X(1,L==L(b)),X(2,L==L(b)),mrkSizeBig,c2,'o','filled','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceAlpha',1);
% plot closest pairs
scatter(X(1,a),X(2,a),250,'s','filled','MarkerEdgeColor',c1,'MarkerFaceColor','w','MarkerFaceAlpha',1,'LineWidth',1.8);
scatter(X(1,b),X(2,b),250,'s','filled','MarkerEdgeColor',c2,'MarkerFaceColor','w','MarkerFaceAlpha',1,'LineWidth',1.8);


end

