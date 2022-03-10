%%
clear variables;
close all;
%%
cd(fileparts(which('test_marble.m')));
addpath(genpath('../'))
dataDir = ['../../../Data/GeometricMarble/'];

%% Options
tstart = tic;
o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');
saveDir = ['../../../Results/tslink/pointclouds/geometricMarble/',dt,'/trial',id,'/'];
%%
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%%
modelType = 'planecylindersphere';
kappa = 10;
epsiNfa = 1;
gricParam.lambda1 = 0;
gricParam.lambda2 = 1;
%%
namePc = 'marble_01_gt.ply';
pc = pcread(fullfile(dataDir,namePc));
%pc =pcdenoise(pc);
pc = pcdownsample(pc,'random',0.002);
%pc = pcdownsample(pc,'gridAverage',0.2);%
%pc = pcdownsample(pc,'nonuniformGridSample',800); %30


%pc =pcdenoise(pc);
normals = pcnormals(pc);
pc = pointCloud(pc.Location, 'normal', normals);
%%

X = [pc.Location'; pc.Normal'];
n = size(X,2);
%%
figure;
subplot(1,2,1);
scatter3(X(1,:),X(2,:),X(3,:),'filled');
axis equal
subplot(1,2,2);
scatter(X(1,:),X(2,:));
axis equal

%%

axis equal;
optsSampling.model = 'planecylindersphere';
optsSampling.sampling = 'nearest';
optsSampling.num_neigh = 400;
optsSampling.quantile = 0.5; % 0.3
optsSampling.m = 10*n; % - 30000;
optsSampling.robust = 'off';
optsSampling.voting = 'gauss';
cardmss = 3; %
optsp = optsSampling;
optsp.model = 'plane';
Sp =  computeResi(X(1:3,:),optsp);
optsc = optsSampling;
optsc.model = 'cylinder';
Sc =  computeResi(X,optsc);
optss = optsSampling;
optss.model = 'sphere';
Ss =  computeResi(X,optss);
%%
tp = tic;
epsi = 1e-3; 
[Sc.P] = resiToP(Sc.R,epsi);
% fprintf('Refitting cylinder..');
% Sc = refitHp(Sc,X,epsi, optsc);
% fprintf('. done\n')
Pc = Sc.P;
%[Pc, isMeaningful_c] = cleansePrefMat(Sc.R, Sc.P, epsi ,kappa,3, epsiNfa);

[Sp.P] = resiToP(Sp.R,epsi);
% fprintf('Refitting plane..');
% Sp = refitHp(Sp,X,epsi, optsp);
% fprintf('. done\n')
Pp = Sp.P;
%Pp(:,sum(Pp,1)>5);
%[Pp, isMeaningful_p] = cleansePrefMat(Sp.R, Sp.P, epsi ,kappa,3, epsiNfa);

[Ss.P] = resiToP(Ss.R,epsi);
Ps = Ss.P;
P =[Pc, Pp, Ps];
tpend = toc(tp);

fprintf('Preference matrix completed in %i minutes.\n',tpend/60);
%% clustering
tc = tic;
%fprintf('Tslink ..');
%t2s = tic;
%C2 = tslnk(P);
%C2 = prune_small_clust(C2,cardmss+1);
%t2 = toc(t2s);
%fprintf('. completed\n');
fprintf('Tslink gric ..');
t3s = tic;
gricParam.sigma = epsi;
C3 = tslnkGric(P,X,modelType,gricParam);
C3 = prune_small_clust(C3,cardmss+1);
t3 = toc(t3s);
tcend = toc(tc);
fprintf('. completed in %d minutes. \n',round(tcend/60));



%%
figure;
optsShowClust.syms = '.';
optsShowClust.mrkrSize = 10;
displayClusters(X,C3, optsShowClust);

%% select model type for clusters
n = size(X,2);
L = C3;
L = prune_small_clust(L,3);
d = 2; % dimension of manifold
r = 3; % dimension of ambient space
kp = 4; % number of parameters for plane
kc = 7; % number of parameters for cylinder
ks = 4; % number of parameters for sphere
T =zeros(n,1); % labels encoding model types
Tmodel = zeros(1, max(L));
Mp = nan(4,max(L));
Mc = nan(7,max(L));
for j = 1:max(L)
    isInClust = L==j;
    Xi = X(:, isInClust);
    
    % fit models
    mp = fit_plane(Xi(1:3,:));
    [w,pc,radius] = fit_cylinder(Xi(1:3,:),[0,0,1]);
    mc = [w;pc;radius];
    ms = fit_sphere(Xi);
    
    % residuals
    rp = res_plane(Xi,mp);
    rc = res_cylinder(Xi,mc);
    rs = res_cylinder(Xi,mc);
    
    % gric 
    
    rSqrp = rp.^2;
    rSqrc = rc.^2;
    rSqrs = rs.^2;
    [gp, ~,~]  = getGricScore(rSqrp,gricParam.sigma, r, d, kp, gricParam.lambda1, gricParam.lambda2);
    [gc, ~,~]  = getGricScore(rSqrc,gricParam.sigma, r, d, kc, gricParam.lambda1, gricParam.lambda2);
    [gs, ~,~]  = getGricScore(rSqrs,gricParam.sigma, r, d, ks, gricParam.lambda1, gricParam.lambda2);
    scores = [gp,gc,gs];
    [~,ind] = min(scores);
    T(isInClust) = ind;
end

%%
figure;
displayClusters(X,T);


%%
% cmap = brewermap(max(L),'Set3');
% cmap = cmap(randperm(max(L)),:);
% %mrkSize = optsShowClust.mrkrSize;
% mrkSize = 8;
% figure; hold all;
% for j  =2:max(L)
%     if(sum(L==j)<10)
%         continue;
%     end
%     if(Tmodel(j)==1)
%         pcPlane = pointCloud(X(1:3,L==j)');
%         model = pcfitplane(pcPlane,epsi);
%        [~, projPoints]=display3dRectangle(model.Parameters(:),X(1:3, L==j),cmap(j,:));
%        %[~, projPoints]=display3dRectangle(Mp(:,j),X(1:3, L==j),cmap(j,:));
%         %scatter3(projPoints(1,:),projPoints(2,:),projPoints(3,:),mrkSize,cmap(j,:),'filled','MarkerEdgeColor','k');
%         %plot3(projPoints(1,:),projPoints(2,:),projPoints(3,:),'.','MarkerSize',mrkSize,'MarkerFaceColor',cmap(j,:),'MarkerEdgeColor',cmap(j,:));
%         s = scatter3(X(1, L==j),X(2, L==j),X(3, L==j),mrkSize,cmap(j,:),'filled');
%     end
%     if(Tmodel(j)==2)
%         pcCyl = pointCloud(X(1:3,L==j)');
%         model = pcfitcylinder(pcCyl,epsi,[0,0,1]);
%         %pcyl = plot(model);
%         set(pcyl,'FaceAlpha',0.9,'FaceColor',cmap(j,:),'EdgeColor',cmap(j,:),'EdgeAlpha',0.9);
%         modelParam  = convertFiniteCylinderToParam(model.Parameters);
%         projectedPoints = projectPointsOnCylinder(X(1:3,L==j),modelParam);
%         scatter3(projectedPoints(1,:),projectedPoints(2,:),projectedPoints(3,:),mrkSize,cmap(j,:),'filled');
%         %scatter3(X(1, L==j),X(2, L==j),X(3, L==j),mrkSize,cmap(j,:),'filled');
%         
%     end
% end
% axis equal;
% axis off;
%%


%%



