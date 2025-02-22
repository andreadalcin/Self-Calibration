clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')


%%
load("ECCV/Data/Seq005_Clip01_Tracks.mat")

num_objects = max(unique(Data.GtLabel));
num_cameras = size(Data.ySparse,3);
debug = false;

Fs = [];
Fo = [];
do = [];
K = [9.842439e+02 0.000000e+00 6.900000e+02; 0.000000e+00 9.808141e+02 2.331966e+02; 0.000000e+00 0.000000e+00 1.000000e+00];

for i = 1:num_cameras-1
    for j = i+5:min(i+10,num_cameras)
        vis_mask = Data.visibleSparse(:,i) & Data.visibleSparse(:,j);

        for o = 1:2
            obj_mask = Data.GtLabel == o;

            file1 = sprintf("ECCV/Sequences/Seq038_Clip02/0000%02d.png", i);
            file2 = sprintf("ECCV/Sequences/Seq038_Clip02/0000%02d.png", j);

            % normalizedPoints1 = normalise2dpts(Data.ySparse(:,vis_mask & obj_mask,i));
            % normalizedPoints2 = normalise2dpts(Data.ySparse(:,vis_mask & obj_mask,j));

            normalizedPoints1 = Data.ySparse(:,vis_mask & obj_mask,i);
            normalizedPoints2 = Data.ySparse(:,vis_mask & obj_mask,j);

            m1 = normalizedPoints1(1:2,:)';
            m2 = normalizedPoints2(1:2,:)';

            if size(m1,1) >= 8 && size(m2,1) >= 8
                fLMedS = fund_lin(m2',m1',[]);
                fLMedS = fund_nonlin(fLMedS,m2',m1');
                % [fLMedS,~] = estimateFundamentalMatrix(m1,m2,'Method','LMedS');
                do = [do; mean(abs(sampson_fund(fLMedS,normalizedPoints1,normalizedPoints2)))];

                Fs(:,:,i,j,o) = fLMedS;
                Fo(:,:,size(Fo,3)+1) = fLMedS;
            end

            if debug
                figure(1);
                I1 = imread(file1);
                imshow(I1)
                axis on
                hold on
                plot(m1(:,1),m1(:,2),'r+','MarkerSize',10)
                epiLines = epipolarLine(fLMedS',m2(:,:));
                points = lineToBorderPoints(epiLines,size(I1));
                line(points(:,[1,3])',points(:,[2,4])');
                hold off

                figure(2);
                I2 = imread(file2);
                imshow(I2)
                axis on
                hold on
                plot(m2(:,1),m2(:,2),'r+','MarkerSize',10)
                epiLines = epipolarLine(fLMedS,m1(:,:));
                points = lineToBorderPoints(epiLines,size(I2));
                line(points(:,[1,3])',points(:,[2,4])');
                hold off
            end
        end
    end
end

do = do ./ sum(do);

Fo(:,:,1) = [];
width = 1.392000e+03;
height = 5.120000e+02;

priors = ones(size(Fo,3),1);
priors = priors ./ sum(priors);


%%
load('Dataset/Fs.mat')
num_cameras = size(Fs,3);

% num_cameras = min(3, num_cameras);
width = 4272;
height = 2848;
K = [];
use_multibody = true;
p = 0.5;

Fo = [];
priors = [];

for i = 1:num_cameras-1
    for j = i+1:num_cameras
        if rand > p
            o1 = 1;
            o2 = 2;
            matches1 = pointMatchesInliers1;
            matches2 = pointMatchesInliers2;
        else
            o1 = 2;
            o2 = 1;
            matches1 = pointMatchesInliers2;
            matches2 = pointMatchesInliers1;
        end

        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o1) / norm(Fs(:,:,i,j,o1));
        priors = [priors; matches1(i,j)];
        if use_multibody
            Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,o2) / norm(Fs(:,:,i,j,o2));
            priors = [priors; matches2(i,j)];
        end
    end
end
Fo(:,:,1) = [];

% priors = priors ./ sum(priors);
priors = ones(size(Fo,3),1);
priors = priors ./ sum(priors);


%% Estimate initial focal length and distribution
% f = [];
% for i = 1:size(Fo,3)
%     f0 = getInitialEstimateOfFocalLength(Fo(:,:,i), width, height);
%     f = [f; f0 .* priors(i) / size(f0,1)];
% end
% 
% % % Filter focal lengths not in the peak bin
% h = histogram(f, 50);
% [~, whichbin] = max(h.Values);
% th_low = h.BinEdges(whichbin);
% th_high = h.BinEdges(whichbin + 1);
% 
% f_pct = f;
% f_pct(or(f < th_low, f > th_high)) = [];
% 
% f0 = sum(f_pct) * size(f,1) / size(f_pct,1);
% fprintf("f0: %f\n", f0)


%% Estimate initial focal length and distribution
f = [];
for i = 1:size(Fo,3)
    f0 = getInitialEstimateOfFocalLength(Fo(:,:,i), width, height);
    f = [f; f0];
end

% % Filter focal lengths not in the peak bin
h = histogram(f, 50);
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);

f_pct = f;
f_pct(or(f < th_low, f > th_high)) = [];

f0 = median(f_pct);
fprintf("f0: %f\n", f0)


%%
weights = ones(size(Fo,3),1);
weights = weights ./ sum(weights);

fxs = [];
fys = [];

min_residual = inf;
fx = inf;
fy = inf;

for i = 1:500
    fprintf("%d\n", i);
    [fx_sample, fy_sample, residual] = evaluate(Fo, weights, f0, width, height);
    fxs = [fxs; fx_sample];
    fys = [fys; fy_sample];

    if residual < min_residual
        min_residual = residual;
        fx = fx_sample;
        fy = fy_sample;
    end
end


%%
figure;
h = histogram(fxs,calcnbins(fxs,'fd'));
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);
fx = fxs;
fx(or(fx < th_low, fx > th_high)) = [];
median(fx)
[~,mux] = robustcov(fxs)

figure;
h = histogram(fys,calcnbins(fys,'fd'));
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);
fy = fys;
fy(or(fy < th_low, fy > th_high)) = [];
median(fy)
[~,muy] = robustcov(fys)



function [fx,fy,residual] = evaluate(Fs, weights, f0, width, height)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mu,sigma] = normfit(f);
% f0 = normrnd(mu,sigma);

% [~, ~, ~, outliers] = robustcov(f);
% f(outliers) = [];
% [mu,sigma] = normfit(f);
% f0 = normrnd(mu,sigma);
% fprintf("Sturm: mu, sigma, f0: %f, %f, %f\n", mu, sigma, f0)

K0 = [f0, 0, width / 2; 0, f0, height / 2; 0, 0, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options = optimoptions('lsqnonlin','Display','off', ...
    'Algorithm','levenberg-marquardt',...
    'StepTolerance',1e-20,...
    'FunctionTolerance',1e-20,...
    'MaxIterations',1e2,...
    'MaxFunctionEvaluations',1e6,...
    'TolFun', 1e-20,...
    'TolX',1e-20);

X0 = [K0(1,1) K0(2,2)];

sample = datasample(1:size(Fs,3),3,'Replace',false,'Weights',weights);
subset = Fs(:,:,sample);

loss = @(X) costFunctionMendoncaCipollaFocalOnly(subset, X, width/2, height/2, '2');

K_SK = lsqnonlin(loss, X0, [], [], Options);
K_SK = [K_SK(1) 0 width/2; 0 K_SK(2) height/2; 0 0 1];

fx = K_SK(1,1);
fy = K_SK(2,2);
residual = 0;
% residual = sum(costFunctionMendoncaCipolla(Fs, [K_SK(1,1) K_SK(1,3) K_SK(2,2) K_SK(2,3)], '2'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% d_f0 = (abs(K(1,1) - f0) + abs(K(2,2) - f0)) / 2;
% d_f = (abs(K(1,1) - K_SK(1,1)) + abs(K(2,2) - K_SK(2,2))) / 2;
% d_uv = norm([K(1,3) - K_SK(1,3), K(2,3) - K_SK(2,3)]);
% disp("Number of Fundamental matrices")
% disp(size(Fs,3))
% disp("Initial focal length - error")
% disp(d_f0)
% disp("Focal length - error")
% disp(d_f)
% disp("Principal point - error")
% disp(d_uv)

end
