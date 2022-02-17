clc, clear, close all
format long;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

load("ECCV/Data/Seq038_Clip01_Tracks.mat")

num_objects = max(unique(Data.GtLabel));
num_cameras = size(Data.ySparse,3);

Fs = [];
Fo = [];
K = [9.842439e+02 0.000000e+00 6.900000e+02; 0.000000e+00 9.808141e+02 2.331966e+02; 0.000000e+00 0.000000e+00 1.000000e+00];

for i = 1:num_cameras-1
    for j = i+1:num_cameras
        vis_mask = Data.visibleSparse(:,i) & Data.visibleSparse(:,j);

        for o = 1:1
            obj_mask = Data.GtLabel == o;

            file1 = sprintf("ECCV/Sequences/Seq005_Clip01/0000%02d.png", i);
            file2 = sprintf("ECCV/Sequences/Seq005_Clip01/0000%02d.png", j);

            % normalizedPoints1 = normalise2dpts(Data.ySparse(:,vis_mask & obj_mask,i));
            % normalizedPoints2 = normalise2dpts(Data.ySparse(:,vis_mask & obj_mask,j));

            normalizedPoints1 = Data.ySparse(:,vis_mask & obj_mask,i);
            normalizedPoints2 = Data.ySparse(:,vis_mask & obj_mask,j);

            m1 = normalizedPoints1(1:2,:)';
            m2 = normalizedPoints2(1:2,:)';

            % fLMedS = fund_lin(m2',m1',[]);
            % fLMedS = fund_nonlin(fLMedS,m2',m1');
            [fLMedS,~] = estimateFundamentalMatrix(m1,m2,'Method','LMedS');
            
            Fs(:,:,i,j,o) = fLMedS;
            Fo(:,:,size(Fo,3)+1) = fLMedS;

%             figure;
%             I1 = imread(file1);
%             imshow(I1)
%             axis on
%             hold on
%             plot(m1(:,1),m1(:,2),'r+','MarkerSize',10)
%             % epiLines = epipolarLine(fLMedS',m2(:,:));
%             % points = lineToBorderPoints(epiLines,size(I1));
%             % line(points(:,[1,3])',points(:,[2,4])');
%           hold off

%            figure;
%             I2 = imread(file2);
%             imshow(I2)
%             axis on
%             hold on
%             plot(m2(:,1),m2(:,2),'r+','MarkerSize',10)
%             % epiLines = epipolarLine(fLMedS,m1(:,:));
%             % points = lineToBorderPoints(epiLines,size(I2));
%             % line(points(:,[1,3])',points(:,[2,4])');
%            hold off
        end
    end
end

Fo(:,:,1) = [];

%%
width = 1.392000e+03;
height = 5.120000e+02;
pct0 = 100;
num_bins = 100;

disp("Single body")
evaluate(Fo, [], K, width, height, num_bins, pct0);



function [d_f0, d_f, d_uv] = evaluate(Fs, outliers, K, width, height, num_bins, pct0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [];
for i = 1:size(Fs,3)
    f0 = getInitialEstimateOfFocalLength(Fs(:,:,i), width, height);
    f = [f; f0];
end

figure;
h = histogram(f, num_bins);

% Filter focal lengths not in the peak bin
[~, whichbin] = max(h.Values);
th_low = h.BinEdges(whichbin);
th_high = h.BinEdges(whichbin + 1);
f(or(f < th_low, f > th_high)) = [];

f0 = median(f);
fprintf("f0: %d\n", f0)
K0 = [f0, 0, width / 2; 0, f0, height / 2; 0, 0, 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% residual = [];
% for i = 1:size(Fs,3)
%     f = getInitialEstimateFromF0(Fs(:,:,i), f0, width, height);
%     if isempty(f)
%         residual = [residual; 1e5];
%     else
%         residual = [residual; abs(f0 - f)];
%     end
% end

% Fo = [];
% pct = pct0;
% while (size(Fo,3) <= 6)
%     th_low = prctile(residual,pct);
%     Fo = Fs(:,:,residual <= th_low);
%     pct = pct + 5;
% end
%
% figure;
% hold on
% plot(residual, 'r.');
% yline(th_low);
% hold off

% disp("Outlier removal percentage:")
% o = residual > th_low;
% sum(o(outliers)) / size(outliers,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt','TolFun', 1e-10,'TolX',1e-10,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);

K0 = real([f0, 0, width / 2; 0, f0, height / 2; 0, 0, 1])
K0 = K;
X0 = [K0(1,:) K0(2,2:3)];
K_SK = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fs, X, '3'), X0, [], [], Options);
% K_SK = robustlsqnonlin(Fs, X0, Options);
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];

disp('Intrinsics - Mendonca&Cipolla')
disp(K_SK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d_f0 = (abs(K(1,1) - f0) + abs(K(2,2) - f0)) / 2;
d_f = (abs(K(1,1) - K_SK(1,1)) + abs(K(2,2) - K_SK(2,2))) / 2;
d_uv = norm([K(1,3) - K_SK(1,3), K(2,3) - K_SK(2,3)]);
disp("Number of Fundamental matrices")
disp(size(Fs,3))
disp("Initial focal length - error")
disp(d_f0)
disp("Focal length - error")
disp(d_f)
disp("Principal point - error")
disp(d_uv)

end
