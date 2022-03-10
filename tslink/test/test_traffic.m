% test the gric idea on the traffic dataset
% cfr ProgX
%% change the folder to the one of the m.file
cd(fileparts(which('test_traffic.m')));
addpath(genpath('../'))
%% Options
clean;
tstart = tic;
o_addOutliers = false;
o_addNoise = false;
o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');

%% Test parameters
for modelType = {'plane'}
    close all;
    modelType = modelType{:};
    saveDir = ['results/traffic/',dt,'/',modelType,'/trial',id,'/'];
    sigma = 0.011;
    cardmss = 3; % the max cardinality possible
    nEpsis = 119;
    kappaSigma = linspace(1,60,nEpsis);
    maxIter = 5;
    kappa = 10;
    epsiNfa = 1;
    gricParam.lambda1 = 1;%0.01;
    gricParam.lambda2 = 2;%0.1;
    %% view pairs names
    
    % three motions
    nameDataset{1} ='cars2_06';
    nameDataset{2} ='cars2_07';
    nameDataset{3} ='cars2B';
    nameDataset{4} ='cars3';
    nameDataset{5} ='cars5';
    nameDataset{6} ='cars9';
    nameDataset{7} ='cars10';
    % two motions
    nameDataset{8} ='cars1';
    nameDataset{9} ='cars2_06_g12';
    nameDataset{10} ='cars2_06_g13';
    nameDataset{11} ='cars2_06_g23';
    nameDataset{12} ='cars2_07_g12';
    nameDataset{13} ='cars2_07_g13';
    nameDataset{14} ='cars2_07_g23';
    nameDataset{15} ='cars2';
    nameDataset{16} ='cars2B_g12';
    nameDataset{17} ='cars2B_g13';
    nameDataset{18} ='cars2B_g23';
    nameDataset{19} ='cars3_g12';
    nameDataset{20} ='cars3_g13';
    nameDataset{21} ='cars3_g23';
    nameDataset{22} ='cars4';
    nameDataset{23} ='cars5_g12';
    nameDataset{24} ='cars5_g13';
    nameDataset{25} ='cars5_g23';
    nameDataset{26} ='cars6';
    nameDataset{27} ='cars7';
    nameDataset{28} ='cars8';
    nameDataset{29} ='cars9_g12';
    nameDataset{30} ='cars9_g13';
    nameDataset{31} ='cars9_g23';
    nameDataset{32} ='cars10_g12';
    nameDataset{33} ='cars10_g13';
    nameDataset{34} ='cars10_g23';
    nameDataset{35} ='kanatani1';
    nameDataset{36} ='kanatani2';
    nameDataset{37} ='truck1';
    nameDataset{38} ='truck2';
    
    
    % articulated
    % three motions
    nameDataset{39} ='two_cranes';
    nameDataset{40} = 'articulated';
    
    % two motions
    nameDataset{41} ='two_cranes_g23';
    nameDataset{42} ='two_cranes_g13';
    nameDataset{43} ='two_cranes_g12';
    nameDataset{44} ='people2';
    nameDataset{45} ='people1';
    nameDataset{46} ='kanatani3';
    nameDataset{47} ='head';
    nameDataset{48} ='articulated_g23';
    nameDataset{49} ='articulated_g13';
    nameDataset{50} ='articulated_g12';
    nameDataset{51} ='arm';
    
    
    nameDataset = nameDataset(:);
    %% preallocation
    nSeq = numel(nameDataset);
    C1 = cell(nSeq,nEpsis, maxIter);
    C2 = cell(nSeq, nEpsis, maxIter);
    C3 = cell(nSeq, nEpsis, maxIter);
    me1 = nan(nSeq, nEpsis, maxIter);
    me2 = nan(nSeq, nEpsis, maxIter);
    me3 = nan(nSeq, nEpsis, maxIter);
    t1 = nan(nSeq, nEpsis, maxIter);
    t2 = nan(nSeq, nEpsis, maxIter);
    t3 = nan(nSeq, nEpsis, maxIter);
    for s = 1:nSeq
        fprintf('Testing %s ...\n',nameDataset{s})
        
        %load(nameDataset{s})
        [data,G] = load_data_h155(nameDataset{s},0);
        X = DataKanatani(data,3);
        % remove repeated point
        [G, orderX] = sort(G);
        X = X(:,orderX);
        n = size(X,2);
        epsisVector = kappaSigma.*sigma;
        
        %% iteration
        for iter = 1:maxIter
            %%
            % hp generation of three model type
            optsSampling.model = modelType;
            optsSampling.sampling = 'localized';
            optsSampling.geo = true;
            optsSampling.m = 6*n;
            optsSampling.robust = 'off';
            optsSampling.voting = 'gauss';
            optsSampling.x84 = false;
            S = computeResi(X,optsSampling);
            
            %% threhsolds
            nEpsis = numel(epsisVector);
            for e = 1:nEpsis
                fprintf('iter %i epsi %i \n',iter,e);
                epsi = epsisVector(e);
                S = refitHp(S,X,epsi, optsSampling);
                [S.P] = resiToP(S.R,epsi,optsSampling.x84);
                [P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,3, epsiNfa);
                %% clustering
                t1s = tic;
                C1{s,e,iter} = tlnkg(P);
                C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss);
                t1(s,e,iter) = toc(t1s);
                t2s = tic;
                C2{s,e,iter} = tslnk(P);
                C2{s,e,iter} = prune_small_clust(C2{s,e,iter},cardmss);
                t2(s,e,iter) = toc(t2s);
                t3s = tic;
                gricParam.sigma = epsi;
                C3{s,e,iter} = tslnkGric(P,X,modelType,gricParam);
                C3{s,e,iter} = prune_small_clust(C3{s,e,iter},cardmss);
                t3(s,e,iter) = toc(t3s);
                %% misclassification error
                me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
                me2(s,e,iter) = computeMissError(C2{s,e,iter},G);
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
            sgtitle('Best and worst results')
            for i = 1:numel(testResults)
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
        
    end
    
    
    
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
    %%
    if(o_dumpTestResults)
        save([saveDir,'summary',id], 'tabME','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','optsSampling');
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
            writecell({'ME mean'}, filename,'Sheet',nameSheet,'Range','A1');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B1');
            writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A2');
            writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B2');
            % misclassification error std
            writecell({'ME std'}, filename,'Sheet',nameSheet,'Range','A81');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B81');
            writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A82');
            writematrix([stdMe;mean(stdMe);median(stdMe)],filename,'Sheet',nameSheet,'Range','B82');
            % time
            writecell({'time'}, filename,'Sheet',nameSheet,'Range','A161');
            writematrix(kappaSigma,filename,'Sheet',nameSheet,'Range','B161');
            writecell([nameDataset;{'mean'};{'median'}], filename,'Sheet',nameSheet,'Range','A162');
            writematrix([meanMe;mean(meanMe);median(meanMe)],filename,'Sheet',nameSheet,'Range','B162');
        end
    end
end
%%
tend = toc(tstart);
fprintf('Total time: %i minutes\n',tend/60);
