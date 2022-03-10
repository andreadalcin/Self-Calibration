
function [missrate,mapped_out]=computeMissRatePrinci(group_out,labels_gt)

% Classified points

% Alignment (only over classified points)
mapped_out = bestMap(labels_gt,group_out);
% compute misclassification rate
missrate = sum(labels_gt ~= mapped_out) / numel(labels_gt);


end
