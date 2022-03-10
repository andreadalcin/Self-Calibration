function [C2] = stabilize_homo_segmentation(X,C1,C2,M2)
%STABILIZE_HOMO_SEGMENTATION 
for i = 1:numel(C2)
 
    if(C1(i)>0 && C2(i)==0)
     
        t =  unique(C2(C1==C1(i)));
        t = t(t>0);
        if(~isempty(t))
            r =res_homography(X(:,i),M2(:,t));
            [u,ind] = min(r);
                C2(i) = t(ind);
            
        end
    end
end
end

