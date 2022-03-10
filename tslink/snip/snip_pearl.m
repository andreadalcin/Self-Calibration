% Snippet to test Pearl

% assume that 
% Sf is sampler output for fundamental
% Sh is sampler output for homography
% Sa is sampler output for affine fundamental

% isMeaningful_f indicates meaningful hypothesis for fundamental
% isMeaningful_h indicates meaningful hypothesis for homography
% isMeaningful_a indicates meaningful hypothesis for affine fundamental



Ef = int32([1-Sf.P(:,isMeaningful_f)]');
Eh = int32([1-Sh.P(:,isMeaningful_h)]');
Ea = int32([1-Sa.P(:,isMeaningful_a)]');


E = [Ef;Eh;Ea];

h = GCO_Create(size(X,2),size(E,1));
GCO_SetDataCost(h,E);
lf = ?;
la = ?;
lh = ?;
label_cost = [lf*ones(1,opts1.m),lh*ones(1,opts_p.m),la*ones(1,opts_p.m)];
GCO_SetLabelCost(h,int32(label_cost));
D = exp(-squareform(pdist(X')).^2);
D = D-diag(diag(D));
GCO_SetNeighbors(h,D);
GCO_SetSmoothCost(h,1*ones(size(E,1)) - diag(ones(1,size(E,1))));
% GCO_SetSmoothCost(h,zeros(size(E,1)));
GCO_Expansion(h);
C_pearl = GCO_GetLabeling(h);
C_pearl= prune_small_clust(C_pearl,cardmss);