%% test motion segmentation on the hopkins 155 dataset using affine space
%  (flats) with different dimension
rng(1986)

clear variables;
close all;
clc;
o_dumpTestResults = true;
o_dumpFigure = true;
o_displayFigure = false;
o_project = false;
%% change the folder to the one of the m.file
cd(fileparts(which('test_h155.m')));

addpath(genpath('../'))

%%
dime = 5; % dimension of projected points

nameSliceDataset = 'Traffic3';

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


for modelType = {'mixedFlats'}
    
    modelType = modelType{:};
    switch modelType
        case 'flat2'
            dimFlat = 2;
            cardmss = dimFlat+1;
            thCard = 2;
        case 'flat3'
            dimFlat = 3;
            cardmss = dimFlat+1;
            thCard = 3;
        case 'mixedFlats'
            thCard = 2;
    end
    close all;
    
    % potresti voler cambiare il folder di salvataggio
    saveDir = ['../../../Results/h155/',nameSliceDataset,'/',dt,'/',modelType,'/trial',id,'/'];
    if(~exist(saveDir,'dir'))
        mkdir(saveDir);
    end
    
    sigma = 0.01;
    
    nEpsis = 30;
    kappaSigma =linspace(1,30,nEpsis);
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
    Clust = cell(nSeq, nEpsis, maxIter);
    silhouette1 = zeros(nSeq, nEpsis, maxIter);
    silhouette2 = zeros(nSeq, nEpsis, maxIter);
    me1 = nan(nSeq, nEpsis, maxIter);
    me2 = nan(nSeq, nEpsis, maxIter);
    me3 = nan(nSeq, nEpsis, maxIter);
    t1 = nan(nSeq, nEpsis, maxIter);
    t2 = nan(nSeq, nEpsis, maxIter);
    t3 = nan(nSeq, nEpsis, maxIter);
    
    for s = 1:nSeq
        nameDataset{s} = erase(nameDataset{s}, '_truth.mat');
        fprintf('Testing %s ...\n',nameDataset{s})
        if(strcmp(nameSliceDataset,'Traffic3'))
            nameDataset{s} = nameDataset{s}(2:end-1);
        end
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
            optsSampling.voting = 'binary';
            optsSampling.x84 = false;
            
            if(strcmp(modelType,'mixedFlats'))
                optsSampling.dim_subspace = 2;
                S2 = computeResi(X,optsSampling);
                optsSampling.dim_subspace = 3;
                S3 = computeResi(X,optsSampling);
                %optsSampling.dim_subspace = 4;
                %S4 = computeResi(X,optsSampling);
            else
                optsSampling.dim_subspace = dimFlat;
                S = computeResi(X,optsSampling);
            end
            
            %% threhsolds
            nEpsis = numel(epsisVector);
            for e = 1:nEpsis
                fprintf('sequence %i iter %i epsi %i \n',s, iter,e);
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
                %P = double(S.P>0);
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
                Clust_raw =  C3{s,e,iter};
                C3{s,e,iter} = prune_small_clust(C3{s,e,iter},thCard);
                t3(s,e,iter) = toc(t3s);
                
                %% compute J-silhouette to select the best threhsold
                Clust{s,e,iter} =  C3{s,e,iter};
                
                if(strcmp(modelType,'mixedFlats'))
                    % recover cluster type and compute residuals
                    % determine the type for each cluster
                    clustClass = nan(1,max(Clust{s,e,iter}));
                    resSilhouette = nan(size(X,2),max(Clust{s,e,iter}));
                    M = cell(1,max(Clust{s,e,iter}));
                    for j = 1:max(Clust{s,e,iter})
                        inliers = X(:,Clust{s,e,iter}==j);
                        %%% compute gric for first flat
                        
                        dimFlat2 = 2;
                        if(size(inliers,2)>dimFlat2)
                        k = 2*dimFlat2;
                        d = dimFlat2;
                        r = size(X,1);
                        flat2 = fit_flat(inliers,dimFlat2);
                        resi2 = res_flat(inliers,flat2);
                        rSqr = resi2.^2;
                        [gScore2,~,~] = getGricScore(rSqr,gricParam.sigma,r,d,k,gricParam.lambda1, gricParam.lambda2);
                        else
                            gScore2 = Inf;
                        end
                        %%% compute gric for second flat
                        dimFlat3 = 3;
                        if(size(inliers,2)>dimFlat3)
                        k = 2*dimFlat3;
                        d = dimFlat3;
                        r = size(X,1);
                        flat3 = fit_flat(inliers,dimFlat3);
                        resi3 = res_flat(inliers,flat3);
                        rSqr = resi3.^2;
                        [gScore3,~,~] = getGricScore(rSqr,gricParam.sigma,r,d,k,gricParam.lambda1, gricParam.lambda2);
                        else
                            gScore3 = Inf;
                        end
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
                    
                    [M,res] = recover_flat( X,Clust{s,e,iter} ,dimFlat);
                    resSilhouette = nan(size(X,2),size(M,2));
                    for j = 1:size(M,2)
                        resSilhouette(:,j) = res_flat(X,M(:,j));
                    end
                end
                
                %%% assign each point to its closer model
                %figure;
                %subplot(1,2,1);
                %displayClusters(X,Clust);
                
                for i = 1:size(X,2)
                    if(Clust{s,e,iter}(i)>0)
                        [dist,indClosest] = min(resSilhouette(i,:));
                        if(dist<epsi && indClosest ~=Clust{s,e,iter}(i))
                            Clust{s,e,iter}(i) = indClosest;
                        end
                    end
                end
                %subplot(1,2,2);
                %displayClusters(X,Clust);
                %%
                tol_a = 1e-3;
                score1 = zeros(size(X,2),1);
                score2 = zeros(size(X,2),1);
                nEstiOutlier = sum(Clust{s,e,iter} ==0);
                totPoints = size(X,2);
                percTolOutlier = 0.5;
                if(max(Clust{s,e,iter})>1 && nEstiOutlier< totPoints*percTolOutlier && max(Clust{s,e,iter})<=8)
                    for i = 1:size(X,2)
                        if(Clust{s,e,iter}(i)==0)
                            score1(i) = 0;
                            score2(i) = 0;
                        else
                            a = resSilhouette(i,Clust{s,e,iter}(i));
                            resSilhouette(i,Clust{s,e,iter}(i)) = inf;
                            b = min(resSilhouette(i,:));
                            score1(i) = (b-a)/epsi;
                            score2(i) = b/max(a,tol_a);
                            
                        end
                    end
                end
                silhouette1(s,e,iter) = mean(score1,'omitnan');
                silhouette2(s,e,iter) = mean(score2,'omitnan');
                
                if(isinf(silhouette1(s,e,iter)) ||isinf(silhouette2(s,e,iter)))
                    keyboard;
                end
                
                if(iter == 1)
%                     figure;
%                     subplot(1,2,1);
%                     displayClusters(X,Clust_raw);
%                     title('raw')
%                     subplot(1,2,2);
%                     displayClusters(X,Clust{s,e,iter});
%                     title('removed small clust');
%                     sgtitle(sprintf("epsi %d, sil1: %.2f, sil2: %.2f", e,silhouette1(s,e,iter),silhouette2(s,e,iter) ));
                end
                %% misclassification error
                %me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
                %fprintf('\t me tlinkg: %.2f\n', 100*me1(s,e,iter));
                %me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
               % me3(s,e,iter) = Misclassification(C3{s,e,iter},G);
                me3(s,e,iter) = Misclassification(Clust{s,e,iter},G);
                fprintf('\t me MLink: %.2f\n', 100*me3(s,e,iter));
                fprintf('\t silhouette1: %.2f\n', silhouette1(s,e,iter));
                fprintf('\t silhouette2: %.2f\n', silhouette2(s,e,iter));
            end
                if(o_dumpTestResults)
                    save([saveDir,'summary',id,'_iter',num2str(iter)],'me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','silhouette1','silhouette2','optsSampling');
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
        %%   save silhouette index
        mm = squeeze(mean(100*me3(s,:,:),3));
        ss1 = squeeze(mean(silhouette1(s,:,:),3));
        ss2 = squeeze(mean(silhouette2(s,:,:),3));
        figure;
        plot(epsisVector, mm,'r','LineWidth',2);
        hold on;
        plot(epsisVector,ss1,'b')
        plot(epsisVector,ss2,'g')
        legend('me','(b-a)/\epsilon','b/a')
        silFolder = [saveDir,'Sil/'];
        if(~exist(silFolder,'dir'))
            mkdir(silFolder);
        end
        saveas(gcf,[silFolder,nameDataset{s},id,'_sil'],'epsc');
        fprintf('Finished sequence %d\n',s);
        %%
        close(gcf);
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
    meanSil1 = squeeze(mean(silhouette1,3));
    meanSil2 = squeeze(mean(silhouette2,3));
    meanMeOur = squeeze(mean(me3,3));
    for i = 1:nSeq
        [~,ind1] = max(meanSil1(i,:));
        bestSilMe1(i,1) = meanMeOur(i,ind1);
        [~,ind2] = max(meanSil2(i,:));
        bestSilMe2(i,1) = meanMeOur(i,ind2);
    end
    %%
    if(o_dumpTestResults)
        save([saveDir,'summary',id], 'tabME','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','silhouette1','silhouette2','optsSampling');
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
        writecell({'Silhouette 1 mean'}, filename,'Sheet','Silhouette1','Range','A1');
        writecell(nameDataset, filename,'Sheet','Silhouette1','Range','A2');
        writematrix(meanSil1,filename,'Sheet','Silhouette1','Range','B2');
        
        writecell({'Silhouette 2 mean'}, filename,'Sheet','Silhouette2','Range','A1');
        writecell(nameDataset, filename,'Sheet','Silhouette2','Range','A2');
        writematrix(meanSil2,filename,'Sheet','Silhouette2','Range','B2');
        
        writecell(nameDataset, filename,'Sheet','BestSilhouette','Range','A1');
        writematrix(bestSilMe1,filename, 'Sheet','BestSilhouette','Range','B1');
        writematrix(bestSilMe2,filename, 'Sheet','BestSilhouette','Range','C1');
    end
end
%%
tend = toc(tstart);
fprintf('Total time: %i minutes\n',tend/60);
