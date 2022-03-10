function h = displayPreferredLines(id, X, S, radius)
%DISPLAYPREDERREDLINES
%%
col = [0.5,0.5,0.5]; % color of preferred lines
%%
h = hggroup('Parent', gca);
point = X(:,id);
plot(X(1,id),X(2,id),'ro','MarkerFace','r','MarkerEdge','k','Parent',h);
drawCircle(point(:)',radius,h);

% model
for j = 1:size(S.P,2)
    if(S.P(id,j))
        a = S.H(1,j);
        b = S.H(2,j);
        c = S.H(3,j);
        tol = 1e-3;
        if(abs(b)<tol)
            % vertical line
            slope = Inf;
            intercpt = -c/a;
        else
            slope = -a/b;
            intercpt = -c/b;
        end
        % intersection with circle
        [xout,yout] = linecirc(slope,intercpt,point(1),point(2),radius);
        line(xout,yout,'Color',col,'Linewidth',0.5,'Parent',h);
    
    end
end
plot(X(1,id),X(2,id),'ro','MarkerFace','r','MarkerEdge','k','Parent',h);
axis equal;
end



