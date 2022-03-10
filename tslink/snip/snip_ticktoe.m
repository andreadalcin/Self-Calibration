%% generate data
addpath(genpath('../'))
clear 'variables';
close all;
n = 30;
m  = 50;
u = linspace(0,1,n);
A = [u; u];
B = [u;-u+1];
r = 0.5;
t = linspace(0,2*pi,m);
cx = -0.9;
cy = +0.5;
C = [r*cos(t)+cx;r*sin(t)+cy];

% figure; hold all;
% scatter(A(1,:),A(2,:));
% scatter(B(1,:),B(2,:));
% scatter(C(1,:),C(2,:));
% axis equal;
% axis off;


X = [A,B,C];
X = X +1e-3*randn(size(X));
X= addOutliersInBB(X,20)
gricParam.lambda1 = 0;
gricParam.lambda2 = 2;
gricParam.sigma = 1;
% preference embedding
modelType = 'lc';
optsSampling.model = modelType;
optsSampling.sampling = 'localized';
optsSampling.m = 100; %6
optsSampling.robust = 'off';
optsSampling.voting = 'gauss';
cardmss = 3; % the max cardinality between the one for line and circle
optsl = optsSampling;
optsl.model = 'line';
Sl = computeResi(X,optsl);
optsc = optsSampling;
optsc.model = 'circle';
Sc =  computeResi(X,optsc);

X = X + 1e-2*randn(size(X));
figure; 
scatter(X(1,:),X(2,:));
axis equal;
axis off;

epsi = 1e-2;
gricParam.sigma = epsi;
%Sl = refitHp(Sl,X,epsi, optsl);
[Sl.P] = resiToP(Sl.R,epsi);
%[Pl, isMeaningful_l] = cleansePrefMat(Sl.R, Sl.P, epsi ,kappa,2, epsiNfa);
%Sc = refitHp(Sc,X,epsi, optsc);
[Sc.P] = resiToP(Sc.R,epsi);
%[Pc, isMeaningful_c] = cleansePrefMat(Sc.R, Sc.P, epsi ,kappa,3, epsiNfa);
P =[Sl.P,Sc.P];

[C, Z]= tslnkGric(P,X,modelType,gricParam);


%%
figure;
pltOpts = defaultOptsClustDisplay();
pause(1);
pltOpts.mrkrSize=70;
displayClusters(X,C,pltOpts);


%% plot dendrogram

Z = linkage(double(P>0),'single','jaccard');
figure;
H = dendrogram(Z,0,'orientation','top','ColorThreshold',Z(end-1,3));
set(H,'LineWidth',2)
xticks([]);
set(gca,'FontSize',20);
ylabel('Tanimoto distance')


%// Changing the colours
lineColours = cell2mat(get(H,'Color'));
colourList = unique(lineColours, 'rows');

myColours = flip(brewermap(3,'Set2'));

%// Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines             
for colour = 2:size(colourList,1)
    %// Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    %// Replace the colour for those lines
    lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for line = 1:size(H,1)
    set(H(line), 'Color', lineColours(line,:));
end
%% draw clustering
q = n/5;
G = [ones(n,1); 2*ones(n,1); 3*ones(m,1)];
defaultOpts = defaultOptsClustDisplay();
defaultOpts.syms ='o';
defaultOpts.mrkrSize = 50;
figure;
displayClusters(X,G, defaultOpts);

%% level 1
% small line in 1

N = numel(G);
cmap = brewermap(N,'Set2');
C1 = G;
C1(20:25) = 4;
C1(26:30) = 5;
C1(31:50) = 6;
C1([33,48]) = 6;
C1(2*n+1:2*n+1+m/2) = 7;
% figure;
% displayClusters(X,C1,defaultOpts);
figure; 
displaySLMerge(X, C1,2*n+1,N , 0 ,cmap,[])
saveas(gcf,'merge1','svg');
saveas(gcf,'merge1','fig');
%%


C2 = C1;
%
C2(G==3) = 3;

figure; 
displaySLMerge(X, C2,2*n+1,N , 0 ,cmap,[])
saveas(gcf,'merge2','svg');
saveas(gcf,'merge2','fig');
%%
% figure;
% displayClusters(X,C2,defaultOpts);


figure; 
displaySLMerge(X, C2,1,n , 0 ,cmap,[])
saveas(gcf,'merge3','svg');
saveas(gcf,'merge3','fig');
%%
C3 = C2;
C3(G==1)= 1;
figure; 
displaySLMerge(X, C3,1,n , 0 ,cmap,[])
saveas(gcf,'merge4','svg');
saveas(gcf,'merge4','fig');
%%
%
% figure;
% displayClusters(X,C3,defaultOpts);
figure; 
displaySLMerge(X, C3,n+1,2*n , 0 ,cmap,[])
saveas(gcf,'merge5','svg');
saveas(gcf,'merge5','fig');

%%
C3(G==2) = 2;
figure; 
displaySLMerge(X, C3,n+1,2*n , 0 ,cmap,[])
saveas(gcf,'merge6','svg');
saveas(gcf,'merge6','fig');
%%

figure; 
displaySLMerge(X, G,1,2*n , 0 ,cmap,[])
saveas(gcf,'merge7','svg');
saveas(gcf,'merge7','fig');