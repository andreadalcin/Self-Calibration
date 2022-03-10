% snip to estimate a global sigma for traffic sequence
cd(fileparts(which('snip_traffic_sigma.m')));
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
for mt = {'plane'}
    modelType = mt{:};
    saveDir = ['stats/traffic/',dt,'/trial',id,'/'];
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
    nameDataset = nameDataset(:);
    nSeq = numel(nameDataset);
    sigmaVec_mad = [];
    %% load dataset
    for s = 1:nSeq
        fprintf('Computing simga for %s ...\n',nameDataset{s})
        
        %load(nameDataset{s})
        [data,G] = load_data_h155(nameDataset{s},0);
        X = DataKanatani(data,3);
        [G, orderX] = sort(G);
         [~,resi]=recover_plane(X,G);
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