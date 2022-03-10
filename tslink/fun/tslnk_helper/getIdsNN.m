function idsNN = getIdsNN(queryId ,k, D)
%GETIDSNN L get the id of the k nearest neighbour points of the queryId 
% Note that the query Id is not included in the returned ball ids.

[~, id] = sort(D(queryId,:),'ascend');
idsNN = id(1:k);

end

