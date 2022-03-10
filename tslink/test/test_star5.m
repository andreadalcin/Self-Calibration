% test the gric idea on our beloved star5 dataset
%% change the folder to the one of the m.file
cd(fileparts(which('test_star5.m')));
addpath(genpath('../'))
addpath('../../../Data/synthetic2D');
%% Options
clean;
o_addOutliers = true;
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

maxIter = 50;%50;
kappa = 10;
epsiNfa = 1;
gricParam.lambda1 = 1;
gricParam.lambda2 = 2;
%% prepare dataset

load(nameDataset);
[M,res, mads,sigmas] = recover_line( X,G );
sigma = max(sigmas);

if(o_addOutliers)
    num_out = max(G)*sum(G==1);
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
nEpsis = 9;
kappaSigma = linspace(2,8,nEpsis);
epsisVector = kappaSigma.*sigma;
C1 = cell(nSeq,nEpsis, maxIter);
C2 = cell(nSeq, nEpsis, maxIter);
C3 = cell(nSeq, nEpsis, maxIter);
C4 = cell(nSeq, nEpsis, maxIter);
me1 = nan(nSeq, nEpsis, maxIter);
me2 = nan(nSeq, nEpsis, maxIter);
me3 = nan(nSeq, nEpsis, maxIter);
me4 = nan(nSeq, nEpsis, maxIter);
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
        [P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,cardmss, epsiNfa);
        
        
        %% clustering
        t1s = tic;
        C1{s,e,iter} = tlnkg(P);
        t1(s,e,iter) = toc(t1s);
        C1{s,e,iter} = prune_small_clust(C1{s,e,iter},2);
        %C1{s,e,iter} = keep_k_clust(C1{s,e,iter},5);
        C1{s,e,iter} = prune_outlier_cardgap(C1{s,e,iter},2);
        t2s = tic;
        C2{s,e,iter} = tslnk(P,X);
        t2(s,e,iter) = toc(t2s);
        C2{s,e,iter} = prune_small_clust(C2{s,e,iter},2);
        C2{s,e,iter} = prune_outlier_cardgap(C2{s,e,iter},2);
        
        gricParam.sigma = epsi;
        t3s = tic;
        C3{s,e,iter} = tslnkGric(P,X,'line',gricParam);
        t3(s,e,iter) = toc(t3s);
        C3{s,e,iter} = prune_small_clust(C3{s,e,iter},2);
        C3{s,e,iter} = prune_outlier_cardgap(C3{s,e,iter},2);
        
        t4s = tic;
        C4{s,e,iter} = tlnkGric(P,X,'line',gricParam);
        t4(s,e,iter) = toc(t4s);
        C4{s,e,iter} = prune_small_clust(C4{s,e,iter},2);
        C4{s,e,iter} = prune_outlier_cardgap(C4{s,e,iter},2);
        
        %%
        me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
        me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
        me3(s,e,iter) = computeMissError(C3{s,e,iter},G);
        me4(s,e,iter) = computeMissError(C4{s,e,iter},G);
    end
end

%% package result
testResults(1).name = 'Tlnk';
testResults(2).name = 'TSlnk';
testResults(3).name = 'TSlnkGric';
testResults(4).name = 'TlnkGric';
testResults(1).me = me1;
testResults(2).me = me2;
testResults(3).me = me3;
testResults(4).me = me4;
testResults(1).clust = C1;
testResults(2).clust = C2;
testResults(3).clust = C3;
testResults(4).clust = C4;
testResults(1).timing = t1;
testResults(2).timing = t2;
testResults(3).timing = t3;
testResults(4).timing = t4;
for i = 1:4
    testResults(i).model = modelType;
    testResults(i).epsisVec = epsisVector;
    testResults(i).kappaSigma = kappaSigma;
end
[testResults] = getTestResultsStat(testResults);
%% dump result to file
id = datestr(now,'_HHMMSSFFF');
saveDir = [saveDir,id,'/'];
%%
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
%% visualize result
if(o_displayFigure)
    figure;
    displayMeGraph(testResults, kappaSigma);
    if(o_dumpFigure)
        saveas(gcf,[saveDir,nameDataset, id,'_me.fig']);
        saveas(gcf,[saveDir,nameDataset, id,'_me'],'epsc');
        pause(1)
    end
    %% worst result on the best epsilon
    % save results on different figure
    optsShowClust.mrkrSize = 25;
    optsShowClust.scheme = 'Accent';
    optsShowClust.syms= 'os<d^v>';
    figure;
    for i = 1:numel(testResults)
        figure;
        titleFig = ['best', testResults(i).name];
        displayClustersLine(X(1:2,:),testResults(i).bestMinMeanClust,optsShowClust);
        %title(titleFig);
        if(o_dumpFigure)
            saveas(gcf,[saveDir,nameDataset,id,titleFig,'.fig']);
            saveas(gcf,[saveDir,nameDataset,id,titleFig],'epsc');
        end
        pause(1);
        figure;
        titleFig = ['worst', testResults(i).name];
        displayClustersLine(X(1:2,:),testResults(i).worstMinMeanClust,optsShowClust);
        % title(titleFig);
        if(o_dumpFigure)
            saveas(gcf,[saveDir,nameDataset,id,titleFig,'.fig']);
            saveas(gcf,[saveDir,nameDataset,id,titleFig],'epsc');
        end
        pause(1)
    end
    %
    %     figure;
    %     sgtitle('Best and worst results')
    %     for i = 1:numel(testResults)
    %         subplot(2,3,i)
    %         displayClustersLine(X(1:2,:),testResults(i).bestMinMeanClust,optsShowClust);
    %         title(['best ', testResults(i).name]);
    %         subplot(2,3,3+i)
    %         displayClustersLine(X(1:2,:),testResults(i).worstMinMeanClust,optsShowClust);
    %         title(['worst ', testResults(i).name]);
    %     end
    %     if(o_dumpFigure)
    %         saveas(gcf,[saveDir,nameDataset,id,'_cfr.fig']);
    %         pause(1)
    %     end
end

%% time
timing = zeros(4,1);
for i =1:4
    bestEpsi = testResults(i).bestEpsi;
    timing(i) = mean(testResults(i).timing(1,testResults(i).epsisVec == bestEpsi,:));
end

if(o_dumpTestResults)
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    saveName = [saveDir,nameDataset,id];
    save(saveName,'testResults','X','timing');
end

pause(1);
close all;