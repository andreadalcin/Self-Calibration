function [par,resi] = lmedsfithomography(X)
%LMEDSFITHOMOGRAPHY
n = size(X,2);
num_trial = 10000;
cardmss = 4;
R = nan(n,num_trial);
H = nan(9,num_trial);
score = nan(1,num_trial);
for i = 1:num_trial
    % instantiate mss
    flg_degen = true;
    cont_degen = 0;
    max_cont_degen = 30;
    while(flg_degen && cont_degen<=max_cont_degen)
        %instantiate mss via uniform sampling
        inds = randsample(n, cardmss,false);
        mss = X(:,inds);
        if(~validateMSS_homography(X, inds))
            flg_degen = 1 ;
        else
            h_tentative = fit_homography(mss);
            flg_degen = is_homography_degen(X, h_tentative, inds);
        end
        cont_degen = cont_degen+1;
    end
    
    % create model
    if(flg_degen == 1)
        inds = nan(cardmss,1);
        h1 = rand(9,1);
    else
        h1 = h_tentative;
    end
    
    % compute residuals
    R(:,i) = res_homography(X, h1);
    H(:,i) = h1;
    score(i) = median(R(:,i).^2);
end
[t,u] = min(score);

par  = H(:,u);
resi = R(:,u);

end

