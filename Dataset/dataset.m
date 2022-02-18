clc, clear, close all

load ('Dataset/params.mat')

images = imageDatastore('Dataset/raw');

for i = 1:size(images.Files,1)
    I = images.readimage(i);
    J = undistortImage(I,cameraParams);
    % figure; imshowpair(I,J,'montage');
    % title('Original Image (left) vs. Corrected Image (right)');

    imwrite(J, strrep(images.Files{i},'raw','undistorted'))
end