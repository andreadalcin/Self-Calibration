function [] = displayStemCard(cardVec,numClust)
%DISPLAYSTEMCARD
if(numClust<=numel(cardVec))
h = stem(1:numClust,cardVec(1:numClust));
set(h(1),'Linewidth',1)
end
if(numel(cardVec)>numClust)
    hold on;
    h(2) = stem(numClust+1:numel(cardVec),cardVec(numClust+1:numel(cardVec)));
    set(h(2),'Linewidth',1)
end


end

