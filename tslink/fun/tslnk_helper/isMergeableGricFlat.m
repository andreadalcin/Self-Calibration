function [ok, msScore, msOutput] = isMergeableGricFlat(X, L, i, j, lambda1 , lambda2,sigma, dimFlat)
% Check if two clusters A and B can be merged.
% The test performs the following steps:
% i) a model i on the first cluster is computed
% iia) a model j in the second cluster is computed
% iib) a model ij on the union of the cluster is computed
% iii) the old gric score on i and j is compared against the new girc score of ij
% A cluster is mergeable if the gric after the merge is lower than before
% If clusters can be merged ok = true; otherwise ok = false;
%
%  Input:
%   - X: points vector
%   - L: vector of clusters labels
%   - i and j: index of points belonging to the cluster to be merged
%   -  lambda1, lambda2 gric parameters
%  Output:
%   - ok: tells if the cluster can be merged
%   - gBefore: gric score before the merge
%   - gAfter: gric score after the merge
%   - mi model fitted on the cluster of i
%   - mj model fitted on the cluster of j
%   - mij model fitted on the union of the cluster of i and j

%%------------------------------------------------------------
% gric magic numbers for flat
% cfr Torr for reference
k = 2*dimFlat; % number of parameters (approximated in order to be coherent with plane)
d = dimFlat; % dimension of the manifold
r = size(X,1); % dimenson of the ambient space
%%------------------------------------------------------------

%% precomputations
isInCi = L==L(i);
isInCj = L==L(j);
% consider points in cluster Ci, in cluster Cj and in the union Ci U Cj
Xi = X(:,isInCi);
Xj = X(:,isInCj);
Xij = X(:,isInCi | isInCj);
% fit a model on Ci, Cj and Ci U Cj
mi = fit_flat(Xi,dimFlat);
mj = fit_flat(Xj,dimFlat);
mij = fit_flat(Xij,dimFlat);
% compute the residuals
ri = res_flat(Xi, mi);
rj = res_flat(Xj, mj);
rij = res_flat(Xij, mij);
% compute squared residual
rSqri = ri.^2;
rSqrj = rj.^2;
rSqrij = rij.^2;
if(nargin< 7)
% compute std
    sigmai = std(rSqri);
    sigmaj= std(rSqrj);
    sigmaij = std(rSqrij);
    sigma = min([sigmai, sigmaj, sigmaij]);
end
%% compute gric score
% gric score before the merge (the sum of gric on individual models)
[gi,dfi,mci]  = getGricScore(rSqri,sigma, r, d, k, lambda1, lambda2);
[gj, dfj,mcj] = getGricScore(rSqrj,sigma, r, d, k, lambda1, lambda2);
gBefore = gi+gj;
dfBefore = dfi +dfj;
mcBefore = mci+mcj;
% gric score after the merge
[gAfter,dfAfter,mcAfter]  = getGricScore(rSqrij,sigma, r, d, k, lambda1, lambda2);
%% compare gric score
ok = gAfter < gBefore;
%% package result
msScore.model = strcat('flat',num2str(dimFlat));
msScore.gric.before = gBefore;
msScore.fidelity.before = dfBefore;
msScore.complexity.before = mcBefore;
msScore.gric.after = gAfter;
msScore.fidelity.after = dfAfter;
msScore.complexity.after = mcAfter;

msOutput.Xi = Xi;
msOutput.Xj = Xj;
msOutput.Xij = Xij;
msOutput.mi = mi;
msOutput.mj = mj;
msOutput.mij = mij;
msOutput.ri = ri;
msOutput.rj = rj;
msOutput.rij = rij;
msOutput.isInCi = isInCi;
msOutput.isInCj = isInCj;
end