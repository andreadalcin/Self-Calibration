% test the gric idea on our beloved star5 dataset
%% change the folder to the one of the m.file
cd(fileparts(which('test_circle4.m')));
addpath(genpath('../'))
%% Options
clean;
o_addOutliers = true;
o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
saveDir = ['../../../Results/tslink_test/jdata/circle4/',dt,'/'];
%%
nameDataset = 'circle4';
modelType = 'circle';
cardmss = 3;

maxIter = 50;
kappa = 10;
epsiNfa = 1;
gricParam.lambda1 = 0;
gricParam.lambda2 = 2;
%% prepare dataset
load(nameDataset);
[M,res, mads,sigmas] = recover_circle( X,G );
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
nEpsis = 7;
kappaSigma = linspace(2,8,nEpsis);
epsisVector = kappaSigma.*sigma;
C1 = cell(nSeq,nEpsis, maxIter);
C2 = cell(nSeq, nEpsis, maxIter);
C3 = cell(nSeq, nEpsis, maxIter);
me1 = nan(nSeq, nEpsis, maxIter);
me2 = nan(nSeq, nEpsis, maxIter);
me3 = nan(nSeq, nEpsis, maxIter);
%%
s = 1;
for iter = 1:maxIter
    % hp generation
    optsSampling.model = modelType;
    optsSampling.sampling = 'localized';
    optsSampling.m = 20*n;
    optsSampling.robust = 'off';
    optsSampling.voting = 'gauss';
    cardmss = 3;
    S = computeResi(X,optsSampling);
    
    
    nEpsis = numel(epsisVector);
    for e = 1:nEpsis
        fprintf('epsi %i iteration %i \n',e,iter);
        epsi = epsisVector(e);
        S = refitHp(S,X,epsi, optsSampling);
        [S.P] = resiToP(S.R,epsi);
        [P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,cardmss, epsiNfa);
        %% clustering
        t1s = tic;
        C1{s,e,iter} = tlnkg(P);
        t1(s,e,iter) = toc(t1s);
        C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss);
        C1{s,e,iter} = prune_outlier_cardgap(C1{s,e,iter},cardmss);
        
        t2s = tic;
        C2{s,e,iter} = tslnk(P);
        C2{s,e,iter} = prune_small_clust(C2{s,e,iter},cardmss);
        t2(s,e,iter) = toc(t2s);
        C2{s,e,iter} = prune_outlier_cardgap(C2{s,e,iter},cardmss);
        
        gricParam.sigma = epsi;
        t3s = tic;
        C3{s,e,iter} = tslnkGric(P,X,modelType,gricParam);
        t3(s,e,iter) = toc(t3s);
        C3{s,e,iter} = prune_small_clust(C3{s,e,iter},cardmss);
        C3{s,e,iter} = prune_outlier_cardgap(C3{s,e,iter},cardmss);
        
        %%
        me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
        me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
        me3(s,e,iter) = computeMissError(C3{s,e,iter},G);
    end
end

%% package result
testResults(1).name = 'Tlnk';
testResults(2).name = 'TSlnk';
testResults(3).name = 'TSlnkG';
testResults(1).me = me1;
testResults(2).me = me2;
testResults(3).me = me3;
testResults(1).clust = C1;
testResults(2).clust = C2;
testResults(3).clust = C3;

testResults(1).timing = t1;
testResults(2).timing = t2;
testResults(3).timing = t3;
for i = 1:3
    testResults(i).model = modelType;
    testResults(i).epsisVec = epsisVector;
end
[testResults] = getTestResultsStat(testResults);
%% dump result to file
id = datestr(now,'_HHMMSSFFF');
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
    end
    %% worst result on the best epsilon
    
    % save results on different figure
    figure;
    for i = 1:numel(testResults)
        figure;
        titleFig = ['best', testResults(i).name];
        displayClustersCircle(X(1:2,:),testResults(i).bestMinMeanClust,optsShowClust);
        %title(titleFig);
        if(o_dumpFigure)
            saveas(gcf,[saveDir,nameDataset,id,titleFig,'.fig']);
            saveas(gcf,[saveDir,nameDataset,id,titleFig],'epsc');
        end
        pause(1);
        figure;
        titleFig = ['worst', testResults(i).name];
        displayClustersCircle(X(1:2,:),testResults(i).worstMinMeanClust,optsShowClust);
        % title(titleFig);
        if(o_dumpFigure)
            saveas(gcf,[saveDir,nameDataset,id,titleFig,'.fig']);
            saveas(gcf,[saveDir,nameDataset,id,titleFig],'epsc');
        end
        pause(1)
    end
    
    
    %     figure;
    %     sgtitle('Best and worst results')
    %     for i = 1:numel(testResults)
    %         subplot(2,3,i)
    %
    %         displayClustersCircle(X(1:2,:),testResults(i).bestMinMeanClust,optsShowClust);
    %
    %         title(['best ', testResults(i).name]);
    %         subplot(2,3,3+i)
    %
    %         displayClustersCircle(X(1:2,:),testResults(i).worstMinMeanClust,optsShowClust);
    %
    %         title(['worst ', testResults(i).name]);
    %     end
    %     if(o_dumpFigure)
    %         saveas(gcf,[saveDir,nameDataset,id,'_cfr.fig']);
    %
    %     end
end

%% time
b = zeros(3,1);
for i =1:3
    bestEpsi = testResults(i).bestEpsi;
    b(i) = mean(testResults(i).timing(1,testResults(i).epsisVec == bestEpsi,:));
end

if(o_dumpTestResults)
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    saveName = [saveDir,nameDataset,id];
    save(saveName,'testResults','X','b');
end
