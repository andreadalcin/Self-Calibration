% test the pearl algorithm on adelaide dataset
%% change the folder to the one of the m.file

cd(fileparts(which('test_pearl_adelFM.m')));
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
for modelType = {'fundamental','haf','affine_fundamental'}
    close all;
    modelType = modelType{:};
    saveDir = ['results/adelFM/pearl/',dt,'/',modelType,'/trial',id,'/'];
    sigma = 0.005;
    cardmss = 8; % the max cardinality possible
    nEpsis = 15;%23;
    kappaSigma = linspace(1,12,nEpsis);
    maxIter = 1%5;
    kappa = 10;
    epsiNfa = 1;
    gricParam.lambda1 = 1;%0.01;
    gricParam.lambda2 = 2;%0.1;
    %% view pairs names
    
    nameDataset{1}   = 'biscuit';
    nameDataset{2}   = 'biscuitbook';
    nameDataset{3}   = 'biscuitbookbox';
    nameDataset{4}   = 'boardgame';
    nameDataset{5}   = 'book';
    nameDataset{6}   = 'breadcartoychips';
    nameDataset{7}   = 'breadcube';
    nameDataset{8}   = 'breadcubechips';
    nameDataset{9}   = 'breadtoy';
    nameDataset{10}  = 'breadtoycar';
    nameDataset{11}  = 'carchipscube';
    nameDataset{12}  = 'cube';
    nameDataset{13}  = 'cubebreadtoychips';
    nameDataset{14}  = 'cubechips';
    nameDataset{15}  = 'cubetoy';
    nameDataset{16}  = 'dinobooks';
    nameDataset{17}  = 'game';
    nameDataset{18}  = 'gamebiscuit';
    nameDataset{19}  = 'toycubecar';
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
    for s = 2;%1:nSeq
        fprintf('Testing %s ...\n',nameDataset{s})
        
        load(nameDataset{s})
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
                [Pf, isMeaningful_f] = cleansePrefMat(Sf.R, Sf.P, epsi ,kappa,3, epsiNfa);
                P =[Ph,Pa, Pf];
                %% pearl clustering
                % labels costs
                %----------------
                lf = 12000;
                la = 8000;
                lh = 4000;
                %---------------
                if(0)
                    Ef = int32([1-Sf.P(:,isMeaningful_f)]);
                    Eh = int32([1-Sh.P(:,isMeaningful_h)]);
                    Ea = int32([1-Sa.P(:,isMeaningful_a)]);
                    label_cost = [lf*ones(1,sum(isMeaningful_f)),lh*ones(1,sum(isMeaningful_h)),la*ones(1,sum(isMeaningful_a))];
                else
                    Ef = int32([1-Sf.P]);
                    Eh = int32([1-Sh.P]);
                    Ea = int32([1-Sa.P]);
                    label_cost = [lf*ones(1,size(Ef,2)),lh*ones(1,size(Eh,2)),la*ones(1,size(Ea,2))];
                end
                
                E = [Ef,Eh,Ea]';
                
                h = GCO_Create(size(E,2),size(E,1));
                GCO_SetDataCost(h,E);
                GCO_SetLabelCost(h,int32(label_cost));
                D = squareform(pdist(X'));
                gamma = max(D(:))/4;
                w = round(100*exp(-(D.^2)./2*gamma));
                w =sparse(triu( w - diag(diag(w))));
                GCO_SetNeighbors(h,w);
                GCO_Expansion(h);
                C_pearl = GCO_GetLabeling(h);
                C_pearl= prune_small_clust(C_pearl,cardmss);
                GCO_Delete(h)
                %
                %% clustering
                t1s = tic;
                C1{s,e,iter} = C_pearl;
                C1{s,e,iter} = prune_small_clust(C1{s,e,iter},cardmss);
                t1(s,e,iter) = toc(t1s);
                t2s = tic;
                %% misclassification error
                me1(s,e,iter) = computeMissError(C1{s,e,iter},G);
            end
        end
        
        %% package result on individual dataset
        testResults(1).name = 'Pearl';
        testResults(1).me = me1(s,:,:);
        testResults(1).clust = C1(s,:,:);
        testResults(1).timing = t1(s,:,:);
        
        for i = 1
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
            %             figure;
            %             displayMeGraph(testResults, kappaSigma);
            %             if(o_dumpFigure)
            %                 saveas(gcf,[saveDir,nameDataset{s}, id,'_me.fig']);
            %             end
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
    
    [Pearl_ME, best1] =  min(squeeze(mean(me1,3)),[],2);
    tabME = table(nameDataset,Pearl_ME);
    tabME.Properties.Description = 'Mean Misclassification error attained with the best threhsold';
    
    %%
    Pearl_Time =  mean(squeeze(mean(t1,3)),2);
    tabTime= table(nameDataset,Pearl_Time);
    tabTime.Properties.Description = 'Mean Time';
    %%
    if(o_dumpTestResults)
        %save([saveDir,'summary',id], 'tabME','tabTime','me1','me2','me3','t1','t2','t3','C1','C2','C3','kappaSigma','optsa','optsh','optsf');
        save([saveDir,'summary',id], 'tabME','tabTime','me1','C1','kappaSigma','optsa','optsh','optsf');
    end
    %% dump to datasheet
    if(o_dumpTestResults)
        meVec = {me1,};
        nameSheetVec = {'Pearl'};
        filename = [saveDir,'report_',modelType,id,'.xlsx'];
        writetable(tabME,filename,'Sheet','BestME','Range','A1')
        writecell({'mean'}, filename,'Sheet','BestME','Range','A21');
        writematrix([mean(Pearl_ME)],filename, 'Sheet','BestME','Range','B21');
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