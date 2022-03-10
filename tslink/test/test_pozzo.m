%%
cd(fileparts(which('test_pozzo.m')));
addpath(genpath('../'))
%% Options
clean;
tstart = tic;
o_debugPlot   = true;
o_displayFigure = true;
o_dumpFigure  = true;
o_dumpTestResults = true;
optsShowClust = defaultOptsClustDisplay();
dt = datestr(now,'yyyymmdd');
id = datestr(now,'_HHMMSSFFF');
saveDir = ['results/pointclouds/pozzo/',dt,'/trial',id,'/'];
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

%%
modelType = 'planecylinder';
kappa = 10;
epsiNfa = 1;
gricParam.lambda1 = 1;
gricParam.lambda2 = 2;
%%
%load('pozzo_reduced');
%pc = pointCloud(X);
pc = pcread('pozzo_full_cleaned.ply');
%pc = pointCloud(Points');
pc =pcdenoise(pc);
%pc = pcdownsample(pc,'random',0.5);
%pc = pcdownsample(pc,'gridAverage',0.5);
%pc = pcdownsample(pc,'nonuniformGridSample',5 ); %30

pc =pcdenoise(pc);
normals = pcnormals(pc);
pc = pointCloud(pc.Location, 'normal', normals);

X = [pc.Location'; pc.Normal'];
n = size(X,2);
%%
X(3,:) = -X(3,:);
X(2,:) = -X(2,:);
figure;
subplot(1,2,1);
scatter3(X(1,:),X(2,:),X(3,:),'filled');
axis equal
subplot(1,2,2);
scatter(X(1,:),X(2,:));
axis equal

%%

axis equal;
optsSampling.model = 'planecylinder';
optsSampling.sampling = 'uniform';
optsSampling.quantile = 0.25;
optsSampling.m = 30000; % - 30000;
optsSampling.robust = 'off';
optsSampling.voting = 'gauss';
cardmss = 3; %
optsp = optsSampling;
optsp.model = 'plane';
Sp =  computeResi(X(1:3,:),optsp);
optsc = optsSampling;
optsc.model = 'cylinder';
Sc =  computeResi(X,optsc);
%%
tp = tic;
epsi = 1; %0.5
[Sc.P] = resiToP(Sc.R,epsi);
fprintf('Refitting cylinder..');
%Sc = refitHp(Sc,X,epsi, optsc);
fprintf('. done\n')
Pc = Sc.P;
%[Pc, isMeaningful_c] = cleansePrefMat(Sc.R, Sc.P, epsi ,kappa,3, epsiNfa);

[Sp.P] = resiToP(Sp.R,epsi);
fprintf('Refitting plane..');
%Sp = refitHp(Sp,X,epsi, optsp);
fprintf('. done\n')
Pp = Sc.P;
%Pp(:,sum(Pp,1)>5);
%[Pp, isMeaningful_p] = cleansePrefMat(Sp.R, Sp.P, epsi ,kappa,3, epsiNfa);
P =[Pc, Pp];
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
fprintf('. completed in %i minutes. \n',tcend/60);



%%
figure;
optsShowClust.syms = '.';
optsShowClust.mrkrSize = 250;
displayClusters(X,C3, optsShowClust);

%% select model type for clusters
n = size(X,2);
L = C3;
L = prune_small_clust(L,10);
%L = prune_outlier_cardgap(L,3);
d = 2; % dimension of manifold
r = 3; % dimension of ambient space
kp = 4; % number of parameters for cylinder
kc = 7; % number of parameters for plane
T =zeros(n,1); % labels encoding model types
Tmodel = zeros(1, max(L));
Mp = nan(4,max(L));
Mc = nan(7,max(L));
for j = 1:max(L)
    isInClust = L==j;
    Xi = X(:, isInClust);
    mp = fit_plane(Xi(1:3,:));
    [w,pc,radius] = fit_cylinder(Xi(1:3,:));
    mc = [w;pc;radius];
    rp = res_plane(Xi,mp);
    rc = res_cylinder(Xi,mc);
    rSqrp = rp.^2;
    rSqrc = rc.^2;
    [gp, dfp,mcp]  = getGricScore(rSqrp,gricParam.sigma, r, d, kp, gricParam.lambda1, gricParam.lambda2);
    [gc, dfc,mcc]  = getGricScore(rSqrc,gricParam.sigma, r, d, kc, gricParam.lambda1, gricParam.lambda2);
    if(o_debugPlot & false)
        clf;
        figure(101); hold all;
        scatter3(X(1,:),X(2,:),X(3,:));
        plot3(X(1,isInClust),X(2,isInClust),X(3,isInClust),'ro');
        view([3,8]);
        disp('---------------------------------')
        fprintf('gric plane %d cyl %d\n',gp,gc);
        fprintf('df   plane %d cyl %d\n',dfp,dfc);
        fprintf('mc   plane %d cyl %d\n',mcp,mcc);
        disp('---------------------------------')
        pause
    end
    
    
    if(gc<gp)
        Tmodel(j) = 2;
        T(isInClust) =2;
        Mc(:,j) = mc;
    else
        Tmodel(j) = 1;
        T(isInClust)=1;
        Mp(:,j) = mp;
    end
end
%%
cmap = brewermap(max(L),'Set3');
cmap = cmap(randperm(max(L)),:);
%mrkSize = optsShowClust.mrkrSize;
mrkSize = 15;
figure; hold all;
for j  =1:max(L)
    if(Tmodel(j)==1)
        pcPlane = pointCloud(X(1:3,L==j)');
        model = pcfitplane(pcPlane,epsi/2);
        [~, projPoints]=display3dRectangle(model.Parameters(:),X(1:3, L==j),cmap(j,:));
        %[~, projPoints]=display3dRectangle(Mp(:,j),X(1:3, L==j),cmap(j,:));
        %scatter3(projPoints(1,:),projPoints(2,:),projPoints(3,:),mrkSize,cmap(j,:),'filled','MarkerEdgeColor','k');
        plot3(projPoints(1,:),projPoints(2,:),projPoints(3,:),'.','MarkerSize',mrkSize,'MarkerFaceColor',cmap(j,:),'MarkerEdgeColor',cmap(j,:));
        %s = scatter3(X(1, L==j),X(2, L==j),X(3, L==j),mrkSize,cmap(j,:),'filled');
    end
    if(Tmodel(j)==2)
        pcCyl = pointCloud(X(1:3,L==j)');
        model = pcfitcylinder(pcCyl,epsi/4);
        pcyl = plot(model);
        set(pcyl,'FaceAlpha',0.9,'FaceColor',cmap(j,:),'EdgeColor',cmap(j,:),'EdgeAlpha',0.9);
        modelParam  = convertFiniteCylinderToParam(model.Parameters);
        projectedPoints = projectPointsOnCylinder(X(1:3,L==j),modelParam);
        scatter3(projectedPoints(1,:),projectedPoints(2,:),projectedPoints(3,:),mrkSize,cmap(j,:),'filled');
        %scatter3(X(1, L==j),X(2, L==j),X(3, L==j),mrkSize,cmap(j,:),'filled');
    end
end
axis equal;
axis off;
%%
figure;
displayClusters(X,T);
%save([saveDir,'pozzo'],'X','L','T','C3','P','Tmodel','n','epsi','gricParam')
%%



%%
% pc = pcread('duomoalignedcleaned.ply');
%
% Y = pc.Location';
% nY = size(Y,2);
% %%
% CC = zeros(nY,1);
% RR = zeros(nY,max(L));
% MMP = [];
% MMC = [];
% for j =1:max(L)
%     if(all(T(L==j)==1))
%         % fit a plane
%         m = fit_plane(X(1:3,L==j));
%         MMP(:,j) = m;
%         RR(:,j) = res_plane(Y,m);
%     else
%         m = fit_cylinder_pc(X(:,L==j),epsi);
%         MMC(:,j) = m;
%         RR(:,j) = res_cylinder(Y,m);
%     end
% end
% %%
% for i = 1:nY
%   [dsort,ind] = sort(RR(i,:),'ascend');
%   dratio = dsort(1)/dsort(2);
%   if(dsort(1)< epsi && (dratio < 0.7))
%       CC(i) = ind(1);
%   end
% end
% %%
%
% cmap = brewermap(max(L),'Paired');
% cmap = cmap(randperm(max(L)),:);
% %mrkSize = optsShowClust.mrkrSize;
% mrkSize = 20;
% figure; hold all;
% for j  =1:max(L)
%     if(Tmodel(j)==1)
%         %pcPlane = pointCloud(X(1:3,L==j)');
%         %model = pcfitplane(pcPlane,epsi/2);
%        %[~, projPoints]=display3dRectangle(model.Parameters(:),X(1:3, L==j),cmap(j,:));
%        [~, projPoints]=display3dRectangle(Mp(:,j),X(1:3, L==j),cmap(j,:));
%         scatter3(projPoints(1,:),projPoints(2,:),projPoints(3,:),mrkSize,cmap(j,:),'filled');
%     end
%     if(Tmodel(j)==2)
%         pcCyl = pointCloud(X(1:3,L==j)');
%         model = pcfitcylinder(pcCyl,epsi,[0,0,1]);
%         pcyl = plot(model);
%
%         set(pcyl,'FaceAlpha',0.5,'FaceColor',cmap(j,:),'EdgeColor',cmap(j,:),'EdgeAlpha',0.5);
%         scatter3(X(1, L==j),X(2, L==j),X(3, L==j),mrkSize,cmap(j,:),'filled');
%     end
% end
% axis equal;
% axis off;
% %%
% figure;
% displayClusters(Y,CC);
%
%
