% snip to estimate a global sigma for the adelaide dataset
cd(fileparts(which('snip_adelH_sigma.m')));
addpath(genpath('../'))
%%
clean;

o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');
%% main loop
cont = 1;
for mt = {'homography','affinity'}
    modelType = mt{:};
    saveDir = ['stats/adelH/',dt,'/trial',id,'/'];
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
        elseif(strcmp(modelType,'affinity'))
            [~,resi, madsh]=recover_affinity_geo(X,G);
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
        filename = [saveDir,'sigmaAdelH_',id,'.xlsx'];
        writecell(nameDataset, filename,'Sheet',['mad_',modelType],'Range','A2');
        writetable(Tmad,filename,'Sheet',['mad_',modelType],'Range','B1')
        writecell(nameDataset, filename,'Sheet',['s_',modelType],'Range','A2');
        writetable(Ts,filename,'Sheet',['s_',modelType],'Range','B1')
        writecell(nameDataset, filename,'Sheet',['q_',modelType],'Range','A2');
        writetable(Tq,filename,'Sheet',['q_',modelType],'Range','B1')
    end
    
    disp('Sigma computed.');
end