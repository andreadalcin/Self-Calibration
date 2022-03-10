function [a,b,dClosest] = findClosestInClusters(D, isInA,isInB)
% Find the closest pairs (a,b) between two clusters A,B
n = size(D,1);
dClosest = Inf;
for i = 1:n
    if(isInA(i))
        for j = 1:n
            if(isInB(j)) % nota se tenessi le coppie ordinate non dovrei scorrere tutta la matrice!
                if(D(i,j)<dClosest)
                    a = i;
                    b = j;
                    dClosest = D(i,j);
                end
            end
        end
    end
end
end
