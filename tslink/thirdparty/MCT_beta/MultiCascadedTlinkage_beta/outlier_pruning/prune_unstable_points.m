function [ F ] = prune_unstable_points( X,C, recover_model, res_model, cardmss, epsi)
%PRUNE_UNSTABLE_POINTS
n = size(X,2);
niter = 10;
F = C;
guard = 1;
while(guard)
    F_old = F;
    for i = 1:n
        if(F(i)>0)
            % find inliers to the cluster of the point
            inl = (F == F(i));
            inl(i) = 0; % leave one out
            if(sum(inl)>cardmss)
                theta = recover_model(X, double(inl));
                d = res_model(X(:,i), theta);
                if(d>epsi) 
                    F(i) = 0; % mark it as an outlier
                end
            else
                F(inl)=0;
            end
            
        end
    end  
  
   guard = ~(all(F_old == F));
end

F = prune_small_clust(F,cardmss);
 
% figure;
% display_clust(X,C); axis off;
% hold on;
% scatter(X(1,F==0),X(2,F==0),'kx');

end

