function [ miss ,  fnr, acc, gtgamma] = computeMissError( c, gt )
%FUCKFUCKME compute the me, the right way


% Da riguardare e correggere (permutazioni sbagliate)


npoints = numel(c);
ngroups = max(gt);
Permutations = perms(1:ngroups);
%Permutations = zeros(1,max(gt));
%pause(1);

miss = nan(size(Permutations,1),1);
for j=1:size(Permutations,1)
    gtperm = zeros(npoints,1);
    gtperm(gt>0) = Permutations(j,gt(gt>0));
    miss(j) = sum( c~=gtperm );
end

[miss,temp] = min(miss);


% gamma is the permutation which remaps the gt labels in order to reduce the me
gamma = Permutations(temp,:);
gtgamma = zeros(npoints,1);
gtgamma(gt>0) = gamma(gt(gt>0));


miss =  miss./npoints;

cc = sum(c(gtgamma==0)~=0); % false negative (actual outliers erroneously deemed as inliers)
aa = sum(c(gtgamma==0)==0); % true positive (outlier detected as outliers)
dd  = sum(c(gtgamma>0)>0);   % true negatives (inliers correctly detected)

fnr = cc/(aa+cc);
acc = (aa+dd)/npoints;


  %figure
  %subplot(2,2,1); myscatter(X,c,[]); title('estimated');
  %subplot(2,2,2); myscatter(X, gt,[]); title('groundtruth');
  %subplot(2,2,3); imagesc(c);
  %subplot(2,2,4); imagesc(gtgamma)


end

