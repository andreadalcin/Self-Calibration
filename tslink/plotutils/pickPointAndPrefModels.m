function [picked,s] = pickPointAndPrefModels(tree,X,S,radius,picked,s,useCircle)
%INTERACTWITHCLUSTERS
%drawn  = false(n,1);
%s = cell(1,numel(picked));

if(nargin <7)
    useCircle = false;
end

button = 1;
while(button==1)   % read ginputs until a mouse right-button occurs
    [x,y,button] = ginput(1);
    if(button ==3)
        return;
    end
    id = knnsearch(tree,[x,y],'k',1);
    hold all;
    if(picked(id))
        set(s{id},'Visible','off');
        picked(id) = false;
    else
        if(useCircle)
            s{id}  = displayPreferredCircles(id, X, S, radius);
        else
             s{id} = displayPreferredLines(id, X, S, radius);
        end
        picked(id) = true;
    end
    
end

end

