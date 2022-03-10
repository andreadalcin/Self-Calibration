% test the gric idea on our beloved floorplan dataset
%% change the folder to the one of the m.file
cd(fileparts(which('test_floor.m')));
addpath(genpath('../'))
%% Options
clean;
o_addOutliers = true;
o_addNoise = true;
o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
optsShowClust.mrkrSize = 50;
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');
saveDir = ['results/synth/floor/',dt,'/trial',id,'/'];
%% Test parameters
modelType = 'lc';  % fit lines and circles
cardmss = 3;
%---------------------------------
nEpsis = 1;
kappaSigma = 9; %linspace(10,11,nEpsis);
%----------------------------------
maxIter = 1; % 10
kappa = 10; % cleansing of the preference matrix (the lower the more clean)
epsiNfa = 1;
gricParam.lambda1 = 1;
gricParam.lambda2 = 2;
%% preallocation
nSeq = 3;
nameDataset = cell(nSeq,1);
C1 = cell(nSeq, nEpsis, maxIter);
C2 = cell(nSeq, nEpsis, maxIter);
C3 = cell(nSeq, nEpsis, maxIter);
me1 = nan(nSeq, nEpsis, maxIter);
me2 = nan(nSeq, nEpsis, maxIter);
me3 = nan(nSeq, nEpsis, maxIter);
t1 = nan(nSeq, nEpsis, maxIter);
t2 = nan(nSeq, nEpsis, maxIter);
t3 = nan(nSeq, nEpsis, maxIter);
%% i test sui dataset buoni sono 1 2 3, 4 ha il gt con 12 modelli
for s=1:nSeq
    nameDataset{s} = ['synth_plan',num2str(s)];
    fprintf('Testing %s ...\n',nameDataset{s})
    %% prepare dataset
    load([nameDataset{s},'_labels']);
    if(s==4)
        t = pi/8; 
        R = [cos(t) -sin(t); sin(t), cos(t)];
        X = R*X;
    end
    if(o_addNoise)
        sigma_gt = 0.1;
        X = X + sigma_gt*randn(size(X));
        sigma = sigma_gt;
    else
        sigma = 0.1;
    end
    
    if(o_addOutliers)
        num_out = floor(size(X,2)/2);
        deltaBB = 0.05;
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
    epsisVector = kappaSigma.*sigma;
   %% iteration
    for iter = 1:maxIter
        %%
        % hp generation of two model type
        optsSampling.model = modelType;
        optsSampling.sampling = 'localized';
        optsSampling.m = 10 *n;
        optsSampling.robust = 'off';
        optsSampling.voting = 'gauss';
        cardmss = 3; % the max cardinality between the one for line and circle
        optsl = optsSampling;
        optsl.model = 'line';
        Sl = computeResi(X,optsl);
        optsc = optsSampling;
        optsc.model = 'circle';
        Sc =  computeResi(X,optsc);
        
        %% threhsolds
        nEpsis = numel(epsisVector);
        for e = 1:nEpsis
            fprintf('iter %i epsi %i \n',iter,e);
            epsi = epsisVector(e);
            Sl = refitHp(Sl,X,epsi, optsl);
            [Sl.P] = resiToP(Sl.R,epsi);
            [Pl, isMeaningful_l] = cleansePrefMat(Sl.R, Sl.P, epsi ,kappa,2, epsiNfa);
            Sc = refitHp(Sc,X,epsi, optsc);
            [Sc.P] = resiToP(Sc.R,epsi);
            [Pc, isMeaningful_c] = cleansePrefMat(Sc.R, Sc.P, epsi ,kappa,3, epsiNfa);
            P =[Pl,Pc];
            %% clustering
            t1s = tic;
            %C1{s,e,iter} = tlnkg(P);
            %C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss+1);
            t1(s,e,iter) = toc(t1s);
            t2s = tic;
            %C2{s,e,iter} = tslnk(P);
            %C2{s,e,iter} = prune_small_clust(C2{s,e,iter},cardmss+1);
            t2(s,e,iter) = toc(t2s);
            t3s = tic;
            gricParam.sigma = epsi;
            C3{s,e,iter} = tslnkGric(P,X,modelType,gricParam);
            C3{s,e,iter} = prune_small_clust(C3{s,e,iter},cardmss+1);
            %C3bis{s,e,iter} = prune_outlier_cardgap(C3{s,e,iter},cardmss);
            C3{s,e,iter} = prune_outlier(X,C3{s,e,iter},modelType,epsi ,kappa, epsiNfa);
            t3(s,e,iter) = toc(t3s);
            %%  misclassification error
            %me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
            %me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
            me3(s,e,iter) = computeMissError(C3{s,e,iter},G);
        end
    end
    
    %% package result on individual dataset
    testResults(1).name = 'Tlnk';
    testResults(2).name = 'TSlnk';
    testResults(3).name = 'TSlnkG';
    testResults(1).me = me1(s,:,:);
    testResults(2).me = me2(s,:,:);
    testResults(3).me = me3(s,:,:);
    testResults(1).clust = C1(s,:,:);
    testResults(2).clust = C2(s,:,:);
    testResults(3).clust = C3(s,:,:);
    testResults(1).timing = t1(s,:,:);
    testResults(2).timing = t2(s,:,:);
    testResults(3).timing = t3(s,:,:);
    for i = 1:3
        testResults(i).model = modelType;
        testResults(i).epsisVec = epsisVector;
    end
    [testResults] = getTestResultsStat(testResults);
    %% dump result to file
    if(o_dumpTestResults)
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        saveName = [saveDir,nameDataset{s},id];
        save(saveName,'testResults');
    end
    %% visualize result
    if(o_displayFigure)
        figure;
        displayMeGraph(testResults, kappaSigma);
        if(o_dumpFigure)
            saveas(gcf,[saveDir,nameDataset{s}, id,'_me.fig']);
        end
        %% worst result on the best epsilon
        figure;
        %sgtitle('Best and worst results')
        for i = 3%1:numel(testResults)
            %subplot(2,3,i)
            %displayClusters(X(1:2,:),testResults(i).bestMinMeanClust,optsShowClust);
            
            %title(['best ', testResults(i).name]);
            %subplot(2,3,3+i)
            displayClusters(X(1:2,:),testResults(i).worstMinMeanClust,optsShowClust);
            
            %title(['worst ', testResults(i).name]);
        end
        if(o_dumpFigure)
            saveas(gcf,[saveDir,nameDataset{s},id,'_cfr.fig']);
            
        end
    end
end
%% save summary statistics per dataset
% mean me
[Tlnk_ME, best1] =  min(squeeze(mean(me1,3)),[],2);
[TSlnk_ME, best2] =  min(squeeze(mean(me2,3)),[],2);
[TSlnkG_ME, best3] =  min(squeeze(mean(me3,3)),[],2);
tabME = table(nameDataset,Tlnk_ME, TSlnk_ME,TSlnkG_ME);
tabME.Properties.Description = 'Mean Misclassification error attained with the best threhsold';
% std me
Tlnk_STD = nan(nSeq,1);
TSlnk_STD = nan(nSeq,1);
TSlnkG_STD = nan(nSeq,1);
std1 = std(me1,1,3);
std2 = std(me2,1,3);
std3 = std(me3,1,3);
for k = 1:nSeq
    Tlnk_STD(k) = std1(k, best1(k));
    TSlnk_STD(k) = std2(k, best2(k));
    TSlnkG_STD(k) = std3(k, best3(k));
end

tabSTD = table(nameDataset,Tlnk_STD, TSlnk_STD,TSlnkG_STD);
tabSTD.Properties.Description = 'Std attained with the best threhsold';
% time
Tlink_Time =  mean(squeeze(mean(t1,3)),2);
TSlnk_Time =  mean(squeeze(mean(t2,3)),2);
TSlnkG_Time =  mean(squeeze(mean(t3,3)),2);
tabTime= table(nameDataset,Tlink_Time, TSlnk_Time,TSlnkG_Time);
tabTime.Properties.Description = 'Mean Time';
%%
if(o_dumpTestResults)
    save([saveDir,'summary',id], 'tabME','tabSTD','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','nameDataset','optsl','optsc');
end

%%
load handel
sound(y,Fs)


%% dump to datasheet
meVec = {me1,me2,me3};
nameSheetVec = {'Tlnk','TSlnk','TslnkG'};
filename = [saveDir,'report_','floor',id,'.xlsx'];

writetable(tabME,filename,'Sheet','BestME','Range','A1')
writecell({'mean'}, filename,'Sheet','BestME','Range','A21');
writematrix([mean(Tlnk_ME),mean(TSlnk_ME), mean(TSlnkG_ME)],filename, 'Sheet','BestME','Range','B21');
for i = 1:numel(meVec)
    meanMe = squeeze(mean(meVec{i},3));
    stdMe = std(meVec{i},0,3);
    nameSheet = nameSheetVec{i};
    % misclassification error
    writecell({'ME mean'}, filename,'Sheet',nameSheet,'Range','A1');
    writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B1');
    writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A2');
    writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B2');
    % misclassification error std
    writecell({'ME std'}, filename,'Sheet',nameSheet,'Range','A31');
    writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B31');
    writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A32');
    writematrix([stdMe;mean(stdMe);median(stdMe)],filename,'Sheet',nameSheet,'Range','B32');
    % time
    writecell({'time'}, filename,'Sheet',nameSheet,'Range','A61');
    writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B61');
    writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A62');
    writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B62');
end
