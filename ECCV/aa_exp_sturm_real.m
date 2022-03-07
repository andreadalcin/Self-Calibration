clc, clear, close all
format longG;

rng default;

addpath('Synthetic/')
addpath('Synthetic/prior/')
addpath('ECCV/')

% dataset = "fountain";
% dataset = "herzjesu";
% dataset = "castle";
% dataset = "01_amiibo_static";
% dataset = "02_amiibo_motion";
% dataset = "04_amiibo_motion";
% dataset = "06_amiibo_static";
% dataset = "10_amiibo_motion";
dataset = "12_amiibo_motion";
% dataset = "KITTI_Seq005";
% dataset = "KITTI_Seq113";

xaxis_limit = true;

switch dataset
    case "*"
        load('Dataset/Sturm/Seq005_Clip02.mat')
        width = 1392;
        height = 512;
        f_gt = mean([9.842439e+02 6.900000e+02]);
    case "fountain"
        load('Dataset/Sturm/fountain.mat')
        width = 3072;
        height = 2048;
        f_gt = 2759.48;
        Fs = pre_process(Fs);
    case "castle"
        load('Dataset/Sturm/castle.mat')
        width = 3072;
        height = 2048;
        f_gt = 2759.48;
        xaxis_limit = false;
        Fs = pre_process(Fs);
    case "herzjesu"
        load('Dataset/Sturm/herzjesu.mat')
        width = 3072;
        height = 2048;
        f_gt = 2759.48;
        Fs = pre_process(Fs);
    case "01_amiibo_static"
        load('Dataset/Sturm/01_amiibo_static.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process(Fs);
    case "02_amiibo_motion"
        load('Dataset/Sturm/02_amiibo_motion.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process(Fs);
    case "04_amiibo_motion"
        load('Dataset/Sturm/04_amiibo_motion.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process(Fs);
    case "05_amiibo_static"
        load('Dataset/Sturm/05_amiibo_static.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process(Fs);
    case "06_amiibo_static"
        load('Dataset/Sturm/06_amiibo_static.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process(Fs);
    case "10_amiibo_motion"
        load('Dataset/Sturm/10_amiibo_motion.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process(Fs);
    case "12_amiibo_motion"
        load('Dataset/Sturm/12_amiibo_motion.mat')
        load('Dataset/params.mat')
        width = cameraParams.Intrinsics.ImageSize(2);
        height = cameraParams.Intrinsics.ImageSize(1);
        f_gt = mean(cameraParams.Intrinsics.FocalLength);
        Fs = pre_process_multi(Fs);
    case "KITTI_Seq005"
        load('Dataset/Sturm/Seq005_Clip01.mat')
        width = 1392;
        height = 512;
        f_gt = mean([9.842439e+02 6.900000e+02]);
        Fs = pre_process(Fs);
    case "KITTI_Seq113"
        load('Dataset/Sturm/Seq113_Clip01.mat')
        width = 1392;
        height = 512;
        f_gt = mean([9.842439e+02 6.900000e+02]);
end

init_focal_length(Fs, width, height, f_gt, dataset, xaxis_limit);


function [f0] = init_focal_length(Fs, width, height, f_gt, dataset, xaxis_limit)

% Output from methods
[mu_1, ~] = sturm1(Fs, width, height)
[mu_1_1, ~] = sturm1_1(Fs, width, height)
[mu_f0, ~, x, ySix] = sturm2(Fs, width, height);
mu_f0

% Plot KDE
figure('Name', 'Sturm plot results', 'NumberTitle', 'off');
% title(sprintf('Rob. Initialization - Synthetic - %d fundamental matrices', num_f), 'FontSize', 20)
hold on
h = plot((x - f_gt) / f_gt, ySix ./ max(ySix),'k-','LineWidth',3);
set(h,'LineSmoothing','On')
set(gca,'FontSize',26);
xlabel('Relative Error (%)', 'FontSize', 32);
ylabel('Density', 'FontSize', 32);
xticks([-1 0 1]);
yticks([]);
if xaxis_limit
    xlim([-1 1])
end
xline((mu_1_1 - f_gt) / f_gt,'r-','LineWidth',5);
xline((mu_1 - f_gt) / f_gt,'g-','LineWidth',5);
xline((mu_f0 - f_gt) / f_gt,'b-','LineWidth',5);
xline(0,'k--','LineWidth',2)
hold off

end


function [Fo, weights] = pre_process(Fs)
Fo = [];
weights = [];
num_cameras = size(Fs,3);
for i = 1:num_cameras-1
    for j = i+1:num_cameras
        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,1) / norm(Fs(:,:,i,j,1));
    end
end
Fo(:,:,1) = [];
end

function [Fo, weights] = pre_process_multi(Fs)
Fo = [];
weights = [];
num_cameras = size(Fs,3);
for i = 1:num_cameras-1
    for j = i+1:num_cameras
        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,1) / norm(Fs(:,:,i,j,1));
        Fo(:,:,size(Fo,3)+1) = Fs(:,:,i,j,2) / norm(Fs(:,:,i,j,2));
    end
end
Fo(:,:,1) = [];
end

function [Fo, weights] = pre_process_kitti(Fs)
Fo = [];
weights = [];
for i = 1:size(Fs,3)-1
    for j = i+1:size(Fs,3)
        if Fs(:,:,i,j) ~= zeros(3,3)
            Fo(:,:,size(Fo,3) + 1) = Fs(:,:,i,j);
        end
    end
end
Fo(:,:,1) = [];
end