clc, clear, close all

load('Synthetic/Data/camera_canon.mat')

matches = zeros(100,4);


for i = 1:num_cams-1
    for j = i+1:num_cams
        for k = 1:motions
            if (rand > outlier_ratio)
                K_curr = K;
            else
                r = (max_outlier_f - min_outlier_f) .* rand(1) + min_outlier_f;
                K_curr = [r 0 width / 2; 0 r height / 2; 0 0 1];
            end

            Fs(:,:,size(Fs,3) + 1) = sampleFundamental(K_curr, sigma);
        end
    end
end