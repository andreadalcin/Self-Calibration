function idsBall = getIdsBall(queryId ,radius, D)
%GETIDBALL get the id of the points inside a ball centered in queryId of
%given radius. Note that the query Id is not included in the returned ball ids.
idx = 1:size(D,1);
idsBall = idx((D(queryId,:)<=radius));

end

