%% test motion segmentation on the hopkins 155 dataset using affine space
%  (flats) with different dimension

clear variables;
close all;
clc;
o_dumpTestResults = true;
o_dumpFigure = true;
o_displayFigure = true;
o_project = false;
%% change the folder to the one of the m.file
cd(fileparts(which('test_h155.m')));
addpath('../../moseg/moSegXuLi/Tools/Support')
addpath(genpath('../'))

%%
dime = 5; % dimension of projected points

nameSliceDataset = 'CB3';

%% Options
tstart = tic;
do_addOutliers = false;
do_addNoise = false;
do_debugPlot   = true;
do_displayFigure = true;
do_dumpFigure  = true;
do_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');

%% Test parameters


for modelType = {'flat2'}%{'flat2','flat3','mixedFlats'}
    
    modelType = modelType{:};
    switch modelType
        case 'flat2'
            dimFlat = 2;
            cardmss = dimFlat+1;
        case 'flat3'
            dimFlat = 3;
            cardmss = dimFlat+1;
        case 'flats'
    end
    close all;
    
    % potresti voler cambiare il folder di salvataggio
    saveDir = ['../../../Results/h155/',nameSliceDataset,'/',dt,'/',modelType,'/trial',id,'/'];
    if(~exist(saveDir,'dir'))
        mkdir(saveDir);
    end
    
    sigma = 0.01; %0.011;
    
    nEpsis = 3;
    kappaSigma = linspace(1,30,nEpsis);
    maxIter = 5;
    kappa = 10;
    epsiNfa = 1;
    gricParam.lambda1 = 1;
    gricParam.lambda2 = 2;
    %%%
    %% view pairs names
    
    data_folder = '../../../Data/h155/';
    addpath(genpath(data_folder));
    
    temp = load([data_folder,'nameList/name',nameSliceDataset,'.mat'],['name',nameSliceDataset]);
    nameDataset = temp.(['name',nameSliceDataset]);
    nameDataset = nameDataset(:);
    
    %% preallocation
    nSeq = numel(nameDataset);
    C1 = cell(nSeq,nEpsis, maxIter);
    C2 = cell(nSeq, nEpsis, maxIter);
    C3 = cell(nSeq, nEpsis, maxIter);
    silhouette = nan(nSeq, nEpsis, maxIter);
    me1 = nan(nSeq, nEpsis, maxIter);
    me2 = nan(nSeq, nEpsis, maxIter);
    me3 = nan(nSeq, nEpsis, maxIter);
    t1 = nan(nSeq, nEpsis, maxIter);
    t2 = nan(nSeq, nEpsis, maxIter);
    t3 = nan(nSeq, nEpsis, maxIter);
    
    for s = 1:nSeq
        nameDataset{s} = erase(nameDataset{s}, '_truth.mat');
        fprintf('Testing %s ...\n',nameDataset{s})
        
        %load(nameDataset{s})
        [data,G] = load_data_h155(nameDataset{s},0);
        % http://www.vision.jhu.edu/data/hopkins155/
        if(o_project)
            X = DataKanatani(data,dime);
        else
            X = data;
        end
        
        % remove repeated point
        [G, orderX] = sort(G);
        X = X(:,orderX);
        n = size(X,2);
        epsisVector = kappaSigma.*sigma;
        
        %% iteration
        for iter = 1:maxIter
            %%
            optsSampling.model = 'flat';
            optsSampling.sampling = 'localized';
            optsSampling.num_neigh = 200;
            optsSampling.geo = true;
            optsSampling.m = 6*n;
            optsSampling.robust = 'off';
            optsSampling.voting = 'gauss';
            optsSampling.x84 = false;
            
            if(strcmp(modelType,'mixedFlats'))
                optsSampling.dim_subspace = 2;
                S2 = computeResi(X,optsSampling);
                optsSampling.dim_subspace = 3;
                S3 = computeResi(X,optsSampling);
                optsSampling.dim_subspace = 4;
                S4 = computeResi(X,optsSampling);
                
            else
                
                optsSampling.dim_subspace = dimFlat;
                S = computeResi(X,optsSampling);
            end
            
            %% threhsolds
            nEpsis = numel(epsisVector);
            for e = 1:nEpsis
                fprintf('iter %i epsi %i \n',iter,e);
                epsi = epsisVector(e);
                if(strcmp(modelType,'mixedFlats'))
                    cardmss = 4;
                    optsSampling.m = 3*n;
                    optsSampling.dim_subspace = 2;
                    S2 = refitHp_subspace(S2,X,epsi, optsSampling);
                    [S2.P] = resiToP(S2.R,epsi,optsSampling.x84);
                    optsSampling.dim_subspace = 3;
                    S3 = refitHp_subspace(S3,X,epsi, optsSampling);
                    [S3.P] = resiToP(S3.R,epsi,optsSampling.x84);
                    optsSampling.dim_subspace = 4;
                    S4 = refitHp_subspace(S4,X,epsi, optsSampling);
                    [S4.P] = resiToP(S2.R,epsi,optsSampling.x84);
                    S.R = [S2.R, S3.R, S4.R];
                    S.P =[S2.P,S3.P,S4.P];
                else
                    S = refitHp_subspace(S,X,epsi, optsSampling);
                    [S.P] = resiToP(S.R,epsi,optsSampling.x84);
                end
                
                P=S.P;
                %[P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,3, epsiNfa);
                %% clustering
                t1s = tic;
                %C1{s,e,iter} = tlnkg(P);
                %C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss);
                t1(s,e,iter) = toc(t1s);
                t2s = tic;
                %C2{s,e,iter} = tslnk(P);
                %C2{s,e,iter} = prune_small_clust(C2{s,e,iter},cardmss);
                t2(s,e,iter) = toc(t2s);
                t3s = tic;
                gricParam.sigma = epsi;
                C3{s,e,iter} = tslnkGric_h155(P,X,modelType,gricParam);
                C3{s,e,iter} = prune_small_clust(C3{s,e,iter},cardmss);
                t3(s,e,iter) = toc(t3s);
                
                %% compute J-silhouette to select the best threhsold
                Clust =  C3{s,e,iter};
                [M,res] = recover_flat( X,Clust ,dimFlat);
                resSilhouette = nan(size(X,2),size(M,2));
                for j = 1:size(M,2)
                    resSilhouette(:,j) = res_flat(X,M(:,j));
                end
                score = zeros(size(X,2),1);
                for i = 1:size(X,2)
                    if(Clust(i)==0)
                        score(i) = nan;
                    else
                        a = resSilhouette(i,Clust(i));
                        resSilhouette(i,Clust(i)) = inf;
                        b = min(resSilhouette(i,:));
                        score(i) = (b-a)/epsi;
                    end
                end
                silhouette(s,e,iter) = mean(score,'omitnan');
                
                
                if(o_displayFigure)
                    %Clust =  C3{s,e,iter};
                    %figure;
                    %gscatter(X(1,:),X(2,:),Clust);
                    figure; scatter(X(1,:),X(2,:), 10,score,'filled')
                end
                %% misclassification error
                %me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
                %fprintf('\t me tlinkg: %.2f\n', 100*me1(s,e,iter));
                %me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
                me3(s,e,iter) = Misclassification(C3{s,e,iter},G);
                fprintf('\t me MLink: %.2f\n', 100*me3(s,e,iter));
                fprintf('\t silhouette: %.2f\n', silhouette(s,e,iter));
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
            sgtitle('Best and worst results')
            for i = 3%1:numel(testResults)
                subplot(2,3,i)
                displayClusters(X(1:2,:),testResults(i).bestMinMeanClust,optsShowClust);
                
                title(['best ', testResults(i).name]);
                subplot(2,3,3+i)
                displayClusters(X(1:2,:),testResults(i).worstMinMeanClust,optsShowClust);
                
                title(['worst ', testResults(i).name]);
            end
            if(o_dumpFigure)
                saveas(gcf,[saveDir,nameDataset{s},id,'_cfr.fig']);
            end
        end
        fprintf('Finished sequence %d\n',s);
    end
    %%
    
    
    %% save summary statistics per dataset
    
    [Tlnk_ME, best1] =  min(squeeze(mean(me1,3)),[],2);
    [TSlnk_ME, best2] =  min(squeeze(mean(me2,3)),[],2);
    [TSlnkG_ME, best3] =  min(squeeze(mean(me3,3)),[],2);
    tabME = table(nameDataset,Tlnk_ME, TSlnk_ME,TSlnkG_ME);
    tabME.Properties.Description = 'Mean Misclassification error attained with the best threhsold';
    
    %%
    Tlink_Time =  mean(squeeze(mean(t1,3)),2);
    TSlnk_Time =  mean(squeeze(mean(t2,3)),2);
    TSlnkG_Time =  mean(squeeze(mean(t3,3)),2);
    tabTime= table(nameDataset,Tlink_Time, TSlnk_Time,TSlnkG_Time);
    tabTime.Properties.Description = 'Mean Time';
    %% best results according to silhouette
    meanSil = squeeze(mean(silhouette,3));
    meanMeOur = squeeze(mean(me3,3));
    for i = 1:nSeq
        [~,ind] = max(meanSil(i,:));
        bestSilMe(i,1) = meanMeOur(i,ind);
    end
    %%
    if(o_dumpTestResults)
        save([saveDir,'summary',id], 'tabME','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','silhouette','optsSampling');
    end
    %% dump to datasheet
    if(o_dumpTestResults)
        meVec = {me1,me2,me3};
        nameSheetVec = {'Tlnk','TSlnk','TslnkG'};
        filename = [saveDir,'report_',modelType,id,'.xlsx'];
        writetable(tabME,filename,'Sheet','BestME','Range','A1')
        writecell({'mean'}, filename,'Sheet','BestME','Range','A21');
        writematrix([mean(Tlnk_ME),mean(TSlnk_ME), mean(TSlnkG_ME)],filename, 'Sheet','BestME','Range','B21');
        for i = 1:numel(meVec)
            meanMe = squeeze(mean(meVec{i},3));
            stdMe = std(meVec{i},0,3);
            nameSheet = nameSheetVec{i};
            % misclassification error
            writecell({'ME mean'}, filename,'Sheet',nameSheet,'Range','A2');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B2');
            writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A3');
            writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B3');
            % misclassification error std
            writecell({'ME std'}, filename,'Sheet',nameSheet,'Range','A84');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B84');
            writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A85');
            writematrix([stdMe;mean(stdMe);median(stdMe)],filename,'Sheet',nameSheet,'Range','B85');
            % time
            writecell({'time'}, filename,'Sheet',nameSheet,'Range','A164');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B164');
            writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A165');
            writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B165');
        end
        
        % write silhouette index
        writecell({'Silhouette mean'}, filename,'Sheet','Silhouette','Range','A1');
        writecell(nameDataset, filename,'Sheet','Silhouette','Range','A2');
        writematrix(meanSil,filename,'Sheet','Silhouette','Range','B2');
      
        writecell(nameDataset, filename,'Sheet','BestSilhouette','Range','A1');
        writematrix(bestSilMe,filename, 'Sheet','BestSilhouette','Range','B1');
    end
end
%%
tend = toc(tstart);
fprintf('Total time: %i minutes\n',tend/60);
