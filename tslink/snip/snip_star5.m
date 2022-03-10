% test the gric idea on our beloved star5 dataset
addpath('../Data/synthetic2D');
addpath(genpath('../tslink'));
%% Options
clean;
o_debugPlot   = true;
o_addNoise    = false;
o_addOutliers = true;
% plotting options
optsShowClust = defaultOptsClustDisplay();
%% prepare dataset
load('star5');
sigma = 0.0075;
epsi = 3.5*sigma; 
epsi = 5*sigma; 
%%
if(o_addOutliers)
    num_out =sum(G==1);
    deltaBB = 0.0;
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
    line([minx, minx+2*epsi],[miny,miny],'Linewidth',10);
    title('Ground truth')
end

%%


% hp generation

optsSampling.epsi = epsi;
optsSampling.model = 'line';
optsSampling.sampling ='uniform';
optsSampling.m = 10*n;
optsSampling.robust = 'off';
optsSampling.voting = 'gauss';
cardmss = 2;
S = sampler_homof(X,optsSampling);

%%
kappa = 3;
epsiNfa = 1;
[P, isMeaningful] = cleansePrefMat(S.R, S.P, epsi ,kappa,cardmss, epsiNfa);
figure;
subplot(1,2,1);
imagesc(S.P)
title('original pref matrix');
subplot(1,2,2)
imagesc(P);
title('cleansed pref matrix');

%%

gricParam.lambda1 = 1;
gricParam.lambda2 = 2;
gricParam.sigma = epsi;
model ='line';
tic
C1 = tlnkg(P);
toc
tic
C2 = tslnk(P);
toc
tic
C3 = tslnkGric(P,X, model, gricParam);
toc
tic
C4 = tlnkGric(P,X, model, gricParam);
toc


figure;
subplot(2,2,1);
displayClustersLine(X,C1, epsi);
title('T linkge')
line([minx, minx+2*epsi],[miny,miny],'Linewidth',10);
subplot(2,2,2);
displayClustersLine(X,C2, epsi);
title('Single Linkage')
subplot(2,2,3);
displayClustersLine(X,C3, epsi);
title('Multi T (single linkage)')
subplot(2,2,4);
displayClustersLine(X,C4, epsi);
title('Multi T (j-linkage)')

%%


