function [] = display_cylinder_ls( X, par, col)
%DISPLAY_CYLINDER 

W = par(1:3);
PC = par(4:6);
r = par(7);




% deterimine two end-points
XC = bsxfun(@minus,X,PC);
projX = bsxfun(@times,(XC'*W)',W);
k = sign(W'*projX).*sqrt(sum(projX.*projX));
pa = PC + min(k).*W;
pb = PC + max(k).*W;
% scatter3(PC(1) + projX(1,:),PC(2) +projX(2,:),PC(3) +projX(3,:),'c')
% scatter3(pa(1),pa(2),pa(3),'ms');
% scatter3(pb(1),pb(2),pb(3),'ms');

%%

model = cylinderModel([pa;pb;r]);
plot(model)
h = findobj(gca,'Type','Surface');
set(h,'FaceAlpha',0.4,'FaceColor',col,'EdgeColor',col)
h = findobj(gca,'Type','Line');

end

