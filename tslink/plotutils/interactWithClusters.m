function [] = interactWithClusters(tree,X,L)
%INTERACTWITHCLUSTERS

button = 1;
labels  = unique(L);
drawn = false(1,numel(labels));
s = cell(1,numel(labels));
while (button==1)   % read ginputs until a mouse right-button occurs
    [x,y,button] = ginput(1);
    if(button ==3)
        return;
    end
    p1 = knnsearch(tree,[x,y],'k',1);
    label = L(p1);
    hold all;
    if(drawn(labels == label))
        set(s{label},'Visible','off'); 
        drawn(labels == label) = false;
    else
        s{label} = scatter(X(1,L==label),X(2,L==label),200,'o');
        drawn(labels == label) = true;
    end
    
end

end

