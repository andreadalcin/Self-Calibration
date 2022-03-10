function [d] = res_subspace(X,M)
%RES_SUBSPACE 

num_subspace = size(M,2);
if(num_subspace==1)
    d = sqrt(sum((M.P*X).^2));
else
    n = size(X,2);
    d = nan(n, num_subspace);
    for i = 1:num_subspace
        d(:,i) =  sqrt(sum((M{i}.P * X).^2));
    end
end
end

