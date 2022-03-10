function C2 = cascaded_tlnkg_fund(X, C1,M1, opts)
%CASCADED_TLNK perform T-linkage on each structurer of C1
%the idea is to use prime^cluster_label new clusters.
% X matrix dxn n data points
% C1 nx1 clustering vector (labels) 
% M1 matrix collecting models' parameters by column
% opts sampling parameters
% C2 nx1 clustering vector (labels corresponding to expanded clusters)
%%


prime= primes(100);
if(max(C1)>100)
    prime = primes(max(C1));
    warning('C1 has many clusters ');
end

c = zeros(numel(C1),max(C1)); % clustering matrix nxk (k = numero di cluster)
L = ones(size(C1)); % perchè poi moltiplichiamo
for i = 1:max(C1)
    opts.fund =  M1(:,i);
    SC = sampler_homof(X(:,C1==i),opts);
    if(any(all(SC.P==1)))
        c(C1==i,i) = 1;
    else
        c(C1==i,i)  = tlnkg(SC.P);
    end
    %  1 0 0  | 1
    %  1 0 0  | 1
    %  2 0 0  | 2
    %  0 1 0  | 3 
    %  0 1 0  | 3
    %  0 0 1  | 4
    L = L.*( prime(i).^c(:,i));
end

    % prime = [2,3,5]
    % 2^1 *3^0 *5^0
    % 2^1 *3^0 *5^0
    % 2^2 *3^0 *5^0
    % 2^0 *3^1 *5^0
    % ...
% If you want to search for models in outliers
% decomment the following if block and comment the last line

if(sum(C1==0)>0)
    opts.model = 'homography';
    S_outlier = sampler_homof(X(:,C1==0),opts);
    C_outlier = tlnkg(S_outlier.P);
    % insert outliers labeling in L;
    I = ones(size(L));
    I(C1==0) = C_outlier;
    L = L.*( prime(i+1).^I);
end
C2 = grp2idx(L);
% C2(C1==0) = C1;



end

