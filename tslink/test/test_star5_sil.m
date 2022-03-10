% test the gric idea on our beloved star5 dataset
%% change the folder to the one of the m.file
cd(fileparts(which('test_star5.m')));
addpath(genpath('../'))
addpath('../../../Data/synthetic2D');
%% Options
clear variables;
o_addOutliers = false;
o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
saveDir = ['../../../Results/tslink_test/jdata/star5/',dt,'/'];

%%
nameDataset = 'star5';
modelType = 'line';

maxIter = 5;%50;
kappa = 10;
epsiNfa = 1;
gricParam.lambda1 = 1;
gricParam.lambda2 = 2;
%% prepare dataset

load(nameDataset);
[M,res, mads,sigmas] = recover_line( X,G );
sigma = max(sigmas);

if(o_addOutliers)
    num_out = sum(G==1);
    deltaBB = 0.0;
    X = addOutliersInBB(X,num_out,deltaBB);
    G = [G;zeros(num_out,1)];
else
    X = X(:,G>0);
    G = G(G>0);
end

[G, orderX] = sort(G);
X = X(:,orderX);
n = size(X,2);
if(o_debugPlot)
    figure;
    displayClusters(X,G,optsShowClust);
    
    minx= min(X(1,:));
    miny = min(X(2,:));
    title('Ground truth')
end

%% preallocation
nSeq = 1;
nEpsis = 20;
kappaSigma = linspace(1,20,nEpsis);
epsisVector = logspace(-2,2,nEpsis);
nEpsis = numel(epsisVector);
C1 = cell(nSeq,nEpsis, maxIter);
C2 = cell(nSeq, nEpsis, maxIter);
C3 = cell(nSeq, nEpsis, maxIter);
C4 = cell(nSeq, nEpsis, maxIter);
me1 = nan(nSeq, nEpsis, maxIter);
me2 = nan(nSeq, nEpsis, maxIter);
me3 = nan(nSeq, nEpsis, maxIter);
me4 = nan(nSeq, nEpsis, maxIter);
sil1 = nan(nSeq, nEpsis, maxIter);
sil2 = nan(nSeq, nEpsis, maxIter);
%%
s = 1;
for iter = 1:maxIter
    % hp generation
    optsSampling.model = modelType;
    optsSampling.sampling = 'localized';
    optsSampling.m = 10*n;
    optsSampling.robust = 'off';
    optsSampling.voting = 'gauss';
    optsSampling.quantile = 0.25;
    cardmss = 2;
    S = computeResi(X,optsSampling);
    
    
    nEpsis = numel(epsisVector);
    for e = 1:nEpsis
        fprintf('epi %i iteration %i \n',e,iter);
        epsi = epsisVector(e);
        S = refitHp(S,X,epsi, optsSampling);
        [S.P] = resiToP(S.R,epsi);
        %[P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,cardmss, epsiNfa);
        P = S.P;
        
        %% clustering
        %         t1s = tic;
        %         C1{s,e,iter} = tlnkg(P);
        %         t1(s,e,iter) = toc(t1s);
        %         C1{s,e,iter} = prune_small_clust(C1{s,e,iter},2);
        %         %C1{s,e,iter} = keep_k_clust(C1{s,e,iter},5);
        %         C1{s,e,iter} = prune_outlier_cardgap(C1{s,e,iter},2);
        %         t2s = tic;
        %         C2{s,e,iter} = tslnk(P,X);
        %         t2(s,e,iter) = toc(t2s);
        %         C2{s,e,iter} = prune_small_clust(C2{s,e,iter},2);
        %         C2{s,e,iter} = prune_outlier_cardgap(C2{s,e,iter},2);
        
        gricParam.sigma = epsi;
        t3s = tic;
        C3{s,e,iter} = tslnkGric(P,X,'line',gricParam);
        t3(s,e,iter) = toc(t3s);
        C3{s,e,iter} = prune_small_clust(C3{s,e,iter},2);
        %
        
        %         t4s = tic;
        %         C4{s,e,iter} = tlnkGric(P,X,'line',gricParam);
        %         t4(s,e,iter) = toc(t4s);
        %         C4{s,e,iter} = prune_small_clust(C4{s,e,iter},2);
        %         C4{s,e,iter} = prune_outlier_cardgap(C4{s,e,iter},2);
        
        %%
        %         me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
        %         me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
    
        %         me4(s,e,iter) = computeMissError(C4{s,e,iter},G);
        %% compute silhouette index
        
        
        Clust = C3{s,e,iter};
        [M,res] = recover_line(X,Clust);
        resSilhouette = nan(size(X,2),size(M,2));
        for j = 1:size(M,2)
            resSilhouette(:,j) = res_line(X,M(:,j));
        end
        %         figure;
        %         subplot(1,2,1);
        %         displayClusters(X,Clust);
        %         title('clustering')
        % modify clustering points are assigned to their closest model
        for i = 1:n
            if(Clust(i)>0)
                [dist,indClosest] = min(resSilhouette(i,:));
                if(dist<epsi && indClosest ~=Clust(i))
                    Clust(i) = indClosest;
                end
            end
        end
        %         subplot(1,2,2);
        %         displayClusters(X,Clust);
        %         title('reassigned')
        
        score1 = zeros(size(X,2),1);
        score2 = zeros(size(X,2),1);
        if(max(Clust)>1)
            for i = 1:size(X,2)
                if(Clust(i)==0)
                    score1(i) = nan;
                    score2(i) = nan;
                else
                    a = resSilhouette(i,Clust(i));
                    resSilhouette(i,Clust(i)) = inf;
                    b = min(resSilhouette(i,:));
                    score1(i) = (b-a)/epsi;
                    score2(i) = b/a;
                end
            end
        end
        sil1(s,e,iter) = mean(score1,'omitnan');
        sil2(s,e,iter) = mean(score2,'omitnan');
        
        C3{s,e,iter} = prune_outlier_cardgap(C3{s,e,iter},2);
        me3(s,e,iter) = computeMissError(Clust,G);
        
    end
end
%%

mm = squeeze(mean(100*me3(s,:,:),3));
ss1 = squeeze(mean(sil1(s,:,:),3));
ss2 = squeeze(mean(sil2(s,:,:),3));
figure;
plot( mm,'r','LineWidth',2);
hold on;
plot(ss1,'b','LineWidth',2)
plot(ss2,'g','LineWidth',2)
legend('me','(b-a)/\epsilon','b/a')
%%

[~,ind1] = max(ss1);
[~,ind2] = max(ss2);
figure;
subplot(1,2,1)
displayClustersLine(X,C3{1,ind1,1},[], epsisVector(ind1));
title('Best silhouette (b-a)/\epsilon')
subplot(1,2,2)
displayClustersLine(X,C3{1,ind2,1},[], epsisVector(ind2));
title('Best silhouette b/a')