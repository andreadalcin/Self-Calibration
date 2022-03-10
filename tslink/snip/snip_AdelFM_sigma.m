% snip to estimate a global sigma for the adelaide dataset
cd(fileparts(which('snip_adelFM_sigma.m')));
addpath(genpath('../'))
%%
clean;

o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = false;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');
%% main loop
cont = 1;
for mt = {'fundamental','homography','affine_fundamental'}
    modelType = mt{:};
    saveDir = ['stats/adelFM/',dt,'/trial',id,'/'];
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
    nSeq = numel(nameDataset);
    sigmaVec_mad = [];
    %% load dataset
    for s = 1:nSeq
        fprintf('Computing simga for %s ...\n',nameDataset{s})
        
        load(nameDataset{s})
        % trow away outliers
        X  = X(:,G>0);
        Y = Y(:,G>0);
        G = G(G>0);
        [G, orderX] = sort(G);
        if(strcmp(modelType,'fundamental'))
            [~,resi, madsf]=recover_fundamental_geo(X,G);
        elseif(strcmp(modelType,'affine_fundamental'))
            [~,resi, madsa]=recover_affine_fundamental_geo(X,G);
        elseif(strcmp(modelType,'homography'))
            [~,resi, madsh]=recover_homography_geo(X,G);
        end
        
        for i = 1:max(G)
            sigmaVec_mad = [sigmaVec_mad, robstd(resi(G==i),'MAD')];
            sigma_mad{s,i} = robstd(resi(G==i),'MAD');
            sigma_s{s,i} = robstd(resi(G==i),'S');
            sigma_q{s,i} = robstd(resi(G==i),'Q');
        end
    end
    
    Tmad = cell2table(sigma_mad);
    Ts = cell2table(sigma_s);
    Tq= cell2table(sigma_q);
    %%
    if(o_displayFigure)
        figure(cont)
        histogram(sigmaVec_mad,10);
       title(['Sigma for ',modelType]);
       cont = cont+1;
    end
    
    
    %% dump to datasheet
    if(o_dumpTestResults)
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        filename = [saveDir,'sigmaAdelFM_',id,'.xlsx'];
        writecell(nameDataset, filename,'Sheet',['mad_',modelType],'Range','A2');
        writetable(Tmad,filename,'Sheet',['mad_',modelType],'Range','B1')
        writecell(nameDataset, filename,'Sheet',['s_',modelType],'Range','A2');
        writetable(Ts,filename,'Sheet',['s_',modelType],'Range','B1')
        writecell(nameDataset, filename,'Sheet',['q_',modelType],'Range','A2');
        writetable(Tq,filename,'Sheet',['q_',modelType],'Range','B1')
    end
    
    disp('Sigma computed.');
end