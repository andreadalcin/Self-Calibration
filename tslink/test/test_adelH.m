% test the gric idea on the line circle and parabola datasets
%% change the folder to the one of the m.file
cd(fileparts(which('test_adelH.m')));
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
for modelType = {'homography'}
    close all;
    modelType = modelType{:};
    saveDir = ['results/adelH/',dt,'/',modelType,'/trial',id,'/'];
    sigma = 0.014;
    cardmss = 4; % the max cardinality possible
    nEpsis = 23;
    kappaSigma = linspace(1,12,nEpsis);
    maxIter = 5;
    kappa = 10;
    epsiNfa = 1;
    gricParam.lambda1 = 1;%0.01;
    gricParam.lambda2 = 2;%0.1;
    %% view pairs names
    
    nameDataset{1}  = 'barrsmith';
    nameDataset{2}  = 'bonhall';
    nameDataset{3}  = 'bonython';
    nameDataset{4}  = 'elderhalla';
    nameDataset{5}  = 'elderhallb';
    nameDataset{6}  = 'hartley';
    nameDataset{7}  = 'johnsona';
    nameDataset{8}  = 'johnsonb';
    nameDataset{9}  = 'ladysymon';
    nameDataset{10} = 'library';
    nameDataset{11} = 'napiera';
    nameDataset{12} = 'napierb';
    nameDataset{13} = 'neem';
    nameDataset{14} = 'nese';
    nameDataset{15} = 'oldclassicswing';
    nameDataset{16} = 'physics';
    nameDataset{17} = 'sene';
    nameDataset{18} = 'unihouse';
    nameDataset{19} = 'unionhouse';
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
        
        load(nameDataset{s})
        % remove repeated point
        [ ~,flg ] = remove_repeated_points(X);
        X = X(:,flg);
        Y = Y(:,flg);
        G = G(flg);
        %X  = X(:,G>0);
        %Y = Y(:,G>0);
        %G = G(G>0);
        [G, orderX] = sort(G);
        %         [~,resi_h, madsh]=recover_homography_geo(X,G);
        %         [~,resi_a, madsa]=recover_affine_fundamental_geo(X,G);
        %         [~,resi_f, madsf]=recover_fundamental_geo(X,G);
        %         sigma_h = mad(resi_h(G>0),1);
        %         sigma_a = mad(resi_a(G>0),1);
        %         sigma_f = mad(resi_f(G>0),1);
        %         sigma = mean([sigma_h, sigma_a, sigma_f]);
        %         if(s ==4)
        %             sigma = min([sigma_h, sigma_a, sigma_f]);
        %         end
        
        X = X(:,orderX);
        Y = Y(:, orderX);
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
            optsSampling.x84 = true;
            optsh = optsSampling;
            optsa = optsSampling;
            optsf = optsSampling;
            % equally distribute the sampled hypotheses
            if(strcmp(modelType,'fundamental'))
                optsa.m = 0;
                optsh.m = 0;
            elseif(strcmp(modelType,'homography'))
                optsf.m = 0;
                optsa.m = 0;
            elseif(strcmp(modelType,'affine_fundamental'))
                optsf.m = 0;
                optsh.m = 0;
            elseif(strcmp(modelType,'hf'))
                optsa.m = 0;
                optsf.m = floor(optsSampling.m/2);
                optsh.m = floor(optsSampling.m/2);
            elseif(strcmp(modelType,'af'))
                optsh.m = 0;
                optsf.m = floor(optsSampling.m/2);
                optsa.m = floor(optsSampling.m/2);
            elseif(strcmp(modelType,'ah'))
                optsf.m = 0;
                optsh.m = floor(optsSampling.m/2);
                optsa.m = floor(optsSampling.m/2);
            elseif(strcmp(modelType,'haf'))
                optsf.m = floor(optsSampling.m/3);
                optsh.m = floor(optsSampling.m/3);
                optsa.m = floor(optsSampling.m/3);
            end
            
            optsh.model = 'homography';
            Sh = computeResi(X,optsh);
            optsa.model = 'affine_fundamental';
            Sa =  computeResi(X,optsa);
            optsf.model = 'fundamental';
            Sf =  computeResi(X,optsf);
            
            %% threhsolds
            nEpsis = numel(epsisVector);
            for e = 1:nEpsis
                fprintf('iter %i epsi %i \n',iter,e);
                epsi = epsisVector(e);
                Sh = refitHp(Sh,X,epsi, optsh);
                [Sh.P] = resiToP(Sh.R,epsi,optsh.x84);
                [Ph, isMeaningful_h] = cleansePrefMat(Sh.R, Sh.P, epsi ,kappa,2, epsiNfa);
                Sa = refitHp(Sa,X,epsi, optsa);
                [Sa.P] = resiToP(Sa.R,epsi,optsh.x84);
                [Pa, isMeaningful_a] = cleansePrefMat(Sa.R, Sa.P, epsi ,kappa,3, epsiNfa);
                Sf = refitHp(Sf,X,epsi, optsf);
                [Sf.P] = resiToP(Sf.R,epsi,optsh.x84);
                [Pf, isMeaningful_p] = cleansePrefMat(Sf.R, Sf.P, epsi ,kappa,3, epsiNfa);
                P =[Ph,Pa, Pf];
                %% clustering
                t1s = tic;
                C1{s,e,iter} = tlnkg(P);
                C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss);
                t1(s,e,iter) = toc(t1s);
                t2s = tic;
                C2{s,e,iter} = tslnk(P);
                C2{s,e,iter} = prune_small_clust(C2{s,e,iter},cardmss);
                t2(s,e,iter) = toc(t1s);
                t3s = tic;
                gricParam.sigma = epsi;
                C3{s,e,iter} = tslnkGric(P,X,modelType,gricParam, img1,Y);
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
                displayImageClusters(Y(1:2,:),testResults(i).bestMinMeanClust,img1,optsShowClust);
                
                title(['best ', testResults(i).name]);
                subplot(2,3,3+i)
                displayImageClusters(Y(1:2,:),testResults(i).worstMinMeanClust,img1,optsShowClust);
                
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
        save([saveDir,'summary',id], 'tabME','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','optsa','optsh','optsf');
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
    end
end
%%
tend = toc(tstart);
fprintf('Total time: %i minutes\n',tend/60);