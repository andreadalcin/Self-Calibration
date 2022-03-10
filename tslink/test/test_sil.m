%% test motion segmentation on the hopkins 155 dataset using affine space
%  (flats) with different dimension
clear variables;
close all;
clc;

datasetName = 'h155'; % adelFM, etc..
nameSliceDataset = 'CB2';


%% change the folder to the one of the m.file
cd(fileparts(which('test_sil.m')));
addpath('../../moseg/moSegXuLi/Tools/Support')
addpath(genpath('../'))

%%
dime = 5; % dimension of projected points

%% Options
tstart = tic;
do_project = false;
do_addOutliers = false;
do_addNoise = false;
do_debugPlot   = true;
do_displayFigure = true;
do_dumpFigure  = true;
do_dumpTestResults = true;
do_cleansePreference = false;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');

%% Test parameters

    
nEpsis = 30;
kappaSigma = linspace(1,30,nEpsis);
sigma = 0.01; 
epsisVector = kappaSigma.*sigma;
maxIter = 5;
kappa = 10;
epsiNfa = 1;
gricParam.lambda1 = 1;
gricParam.lambda2 = 2;

for modelType = {'mixedFlats'}%{'flat2','flat3','mixedFlats'}
    
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
    saveDir = ['../../../Results/',datasetName,'/',nameSliceDataset,'/',dt,'/',modelType,'/trial',id,'/'];
    if(~exist(saveDir,'dir'))
        mkdir(saveDir);
    end
    
    % data folder
    data_folder = ['../../../Data/',datasetName,'/']; % 2DO: come sono caricati i dataset
    addpath(genpath(data_folder));
    
    temp = load([data_folder,'nameList/name',nameSliceDataset,'.mat'],['name',nameSliceDataset]);
    nameSequence = temp.(['name',nameSliceDataset]);
    nameSequence = nameSequence(:);
    
    %% preallocation
    nSeq = numel(nameSequence);
    C1 = cell(nSeq,nEpsis, maxIter);  % T-linkage
    C2 = cell(nSeq, nEpsis, maxIter); % MCT
    C3 = cell(nSeq, nEpsis, maxIter); % Multi-Link
    silhouette1 = nan(nSeq, nEpsis, maxIter);
    silhouette2 = nan(nSeq, nEpsis, maxIter);
    silhouette3 = nan(nSeq, nEpsis, maxIter);
    me1 = nan(nSeq, nEpsis, maxIter);
    me2 = nan(nSeq, nEpsis, maxIter);
    me3 = nan(nSeq, nEpsis, maxIter);
    t1 = nan(nSeq, nEpsis, maxIter);
    t2 = nan(nSeq, nEpsis, maxIter);
    t3 = nan(nSeq, nEpsis, maxIter);

    
    for s = 1:nSeq
        nameSequence{s} = erase(nameSequence{s}, '_truth.mat');
        fprintf('Testing %s ...\n',nameSequence{s})
        
        %load(nameSequence{s})
        [data,G] = load_data_h155(nameSequence{s},0);
        % http://www.vision.jhu.edu/data/hopkins155/
        if(do_project)
            X = DataKanatani(data,dime);
        else
            X = data;
        end
        
        % reorder points only for visualization purpose
        [G, orderX] = sort(G);
        X = X(:,orderX);
        n = size(X,2);
      
        
        %% iteration
        for iter = 1:maxIter
            %% Sampling options
            optsSampling.model = 'flat';
            optsSampling.sampling = 'localized';
            optsSampling.num_neigh = 200;
            optsSampling.geo = true;
            optsSampling.m = 6*n;
            optsSampling.robust = 'off';
            optsSampling.voting = 'gauss';
            optsSampling.x84 = true;
            
            % sample hypotheses
            if(strcmp(modelType,'mixedFlats'))
                optsSampling.dim_subspace = 2;
                S2 = computeResi(X,optsSampling);
                optsSampling.dim_subspace = 3;
                S3 = computeResi(X,optsSampling);
            else
                optsSampling.dim_subspace = dimFlat;
                S = computeResi(X,optsSampling);
            end
            
            %% threhsolds
            nEpsis = numel(epsisVector);
            % define preferences
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
                    %optsSampling.dim_subspace = 4;
                    %S4 = refitHp_subspace(S4,X,epsi, optsSampling);
                    %[S4.P] = resiToP(S2.R,epsi,optsSampling.x84);
                    S.R = [S2.R, S3.R]; %[S2.R, S3.R, S4.R];
                    S.P = [S2.P, S3.P];%[S2.P,S3.P,S4.P];
                else
                    S = refitHp_subspace(S,X,epsi, optsSampling);
                    [S.P] = resiToP(S.R,epsi,optsSampling.x84);
                end
                
                P = S.P;
                if(do_cleansePreference)
                    [P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,3, epsiNfa);
                end

                %% clustering
                % t-linkage
                t1s = tic;
                C1{s,e,iter} = tlnkg(P);
                C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss);
                t1(s,e,iter) = toc(t1s);
                % MCT 
                t2s = tic;
                %C2{s,e,iter} =  MCT da implementare, l'input sono i punti e il clustering C1{s,e,iter} e la soglia
                C2{s,e,iter} = prune_small_clust(C2{s,e,iter},cardmss);
                t2(s,e,iter) = toc(t2s);
                t3s = tic;
                % Multi-link
                gricParam.sigma = epsi;
                C3{s,e,iter} = tslnkGric_h155(P,X,modelType,gricParam);
                C3{s,e,iter} = prune_small_clust(C3{s,e,iter},cardmss);
                t3(s,e,iter) = toc(t3s);
                
                % 2D0: mettere il J-shilouette da qui fino a 249 in una funzione: silhouette1(s,e,iter)= js(X, C1{s,e,iter}, @fit_flat, @res_flat)
                % function s = js(X, Clust, fit_model, res_model), dove fit_model = @fit_flat
                %% compute J-silhouette to select the best threhsold, va fatto per tutti i metodi
                Clust =  C3{s,e,iter};
               
                if(strcmp(modelType,'mixedFlats'))
                    % recover cluster type and compute residuals
                    % determine the type for each cluster
                    clustClass = nan(1,max(Clust));
                    resSilhouette = nan(size(X,2),max(Clust));
                    M = cell(1,max(Clust));
                    for j = 1:max(Clust)
                        inliers = X(:,Clust==j);
                        %%% compute gric for first flat
                        dimFlat2 = 2;
                        k = 2*dimFlat2;
                        d = dimFlat2;
                        r = size(X,1);
                        flat2 = fit_flat(inliers,dimFlat2);
                        resi2 = res_flat(inliers,flat2);
                        rSqr = resi2.^2;
                        [gScore2,~,~] = getGricScore(rSqr,gricParam.sigma,r,d,k,gricParam.lambda1, gricParam.lambda2);
                        
                        %%% compute gric for second flat
                        dimFlat3 = 3;
                        k = 2*dimFlat3;
                        d = dimFlat3;
                        r = size(X,1);
                        flat3 = fit_flat(inliers,dimFlat3);
                        resi3 = res_flat(inliers,flat3);
                        rSqr = resi3.^2;
                        [gScore3,~,~] = getGricScore(rSqr,gricParam.sigma,r,d,k,gricParam.lambda1, gricParam.lambda2);
                        %%% determine model class
                        if(gScore2<gScore3)
                            clustClass(j) = dimFlat2;
                            M{j} = flat2;
                            resSilhouette(:,j) = res_flat(X,flat2);
                        else
                            clustClass(j) = dimFlat3;
                            M{j} = flat3;
                            resSilhouette(:,j) = res_flat(X,flat3);
                        end
                    end
                else  %% compute residual for model of the same class

                    [M,res] = recover_flat( X,Clust ,dimFlat);
                    resSilhouette = nan(size(X,2),size(M,2));
                    for j = 1:size(M,2)
                        resSilhouette(:,j) = res_flat(X,M(:,j));
                    end
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
                silhouette1(s,e,iter) = mean(score,'omitnan');
                
                
                
                if(do_displayFigure)
                    %Clust =  C3{s,e,iter};
                    %figure;
                    %gscatter(X(1,:),X(2,:),Clust);
                    figure; scatter(X(1,:),X(2,:), 10,score,'filled')
                end
                %% misclassification error
                me1(s,e,iter) = Misclassification(C1{s,e,iter},G);
                fprintf('\t me tlinkg: %.2f\n', 100*me1(s,e,iter));
                me2(s,e,iter) = Misclassification(C2{s,e,iter},G);
                fprintf('\t me MCT: %.2f\n', 100*me1(s,e,iter));
                me3(s,e,iter) = Misclassification(C3{s,e,iter},G);
                fprintf('\t me MLink: %.2f\n', 100*me3(s,e,iter));
                fprintf('\t silhouette: %.2f\n', silhouette1(s,e,iter));
            end
        end
        
        %% package result on individual dataset
        testResults(1).name = 'Tlnk';
        testResults(2).name = 'MCT';
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
        if(do_dumpTestResults)
            if ~exist(saveDir, 'dir')
                mkdir(saveDir);
            end
            saveName = [saveDir,nameSequence{s},id];
            save(saveName,'testResults');
        end
        %% visualize result
        if(do_displayFigure)
            figure;
            displayMeGraph(testResults, kappaSigma);
            if(do_dumpFigure)
                saveas(gcf,[saveDir,nameSequence{s}, id,'_me.fig']);
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
            if(do_dumpFigure)
                saveas(gcf,[saveDir,nameSequence{s},id,'_cfr.fig']);
            end
        end
        fprintf('Finished sequence %d\n',s);
    end
    %%
    
    
    %% save summary statistics per dataset
    
    [Tlnk_ME, best1] =  min(squeeze(mean(me1,3)),[],2);
    [TSlnk_ME, best2] =  min(squeeze(mean(me2,3)),[],2);
    [TSlnkG_ME, best3] =  min(squeeze(mean(me3,3)),[],2);
    tabME = table(nameSequence,Tlnk_ME, TSlnk_ME,TSlnkG_ME);
    tabME.Properties.Description = 'Mean Misclassification error attained with the best threhsold';
    
    %%
    Tlink_Time =  mean(squeeze(mean(t1,3)),2);
    TSlnk_Time =  mean(squeeze(mean(t2,3)),2);
    TSlnkG_Time =  mean(squeeze(mean(t3,3)),2);
    tabTime= table(nameSequence,Tlink_Time, TSlnk_Time,TSlnkG_Time);
    tabTime.Properties.Description = 'Mean Time';
    %% best results according to silhouette
    meanSil = squeeze(mean(silhouette,3));
    meanMeOur = squeeze(mean(me3,3));
    for i = 1:nSeq
        [~,ind] = max(meanSil(i,:));
        bestSilMe(i,1) = meanMeOur(i,ind);
    end
    %%
    if(do_dumpTestResults)
        save([saveDir,'summary',id], 'tabME','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','silhouette','optsSampling');
    end
    %% dump to datasheet
    if(do_dumpTestResults)
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
            writecell([nameSequence;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A3');
            writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B3');
            % misclassification error std
            writecell({'ME std'}, filename,'Sheet',nameSheet,'Range','A84');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B84');
            writecell([nameSequence;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A85');
            writematrix([stdMe;mean(stdMe);median(stdMe)],filename,'Sheet',nameSheet,'Range','B85');
            % time
            writecell({'time'}, filename,'Sheet',nameSheet,'Range','A164');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B164');
            writecell([nameSequence;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A165');
            writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B165');
        end
        
        % write silhouette index
        writecell({'Silhouette mean'}, filename,'Sheet','Silhouette','Range','A1');
        writecell(nameSequence, filename,'Sheet','Silhouette','Range','A2');
        writematrix(meanSil,filename,'Sheet','Silhouette','Range','B2');
        
        writecell(nameSequence, filename,'Sheet','BestSilhouette','Range','A1');
        writematrix(bestSilMe,filename, 'Sheet','BestSilhouette','Range','B1');
    end
end
%%
tend = toc(tstart);
fprintf('Total time: %i minutes\n',tend/60);
