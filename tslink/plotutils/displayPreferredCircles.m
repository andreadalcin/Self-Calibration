function h = displayPreferredCircles(id, X, S, radius)
%DISPLAYPREDERREDLINES

col = [0.5,0.5,0.5]; % color of preferred circle
h = hggroup('Parent', gca);
point = X(:,id);
plot(X(1,id),X(2,id),'ro','MarkerFace','r','MarkerEdge','k','Parent',h);
drawCircle(point(:)',radius,h);

% model
for j = 1:size(S.P,2)
    if(S.P(id,j))
        % intersection with circle
        [xout,yout,isInside] = circcirc(point(1),point(2),radius,S.H(1,j),S.H(2,j),S.H(3,j));
        
        
        if(isInside)
            drawCircle(S.H(1:2,j),S.H(3,j),h)
        else

  
            % determine the arc of circumference
            t1 = atan2(yout(1)-S.H(2,j),xout(1)-S.H(1,j));
            t2 = atan2(yout(2)-S.H(2,j),xout(2)-S.H(1,j));
            % map angle between 0,2pi
            if(t1<0)
                t1 = t1+2*pi;
            end
            if(t2<0)
                t2 = t2+2*pi;
            end
            tmax = max(t1,t2);
            tmin = min(t1,t2);
            if(tmax-tmin <=pi)
                thetaStart = tmin;
                thetaEnd = tmax;
            else
                thetaStart = tmax;
                thetaEnd = tmin+2*pi;
            end
            drawArcCircle(S.H(1:2,j),S.H(3,j), thetaStart, thetaEnd, h,col,0.5)
        end
    end
end
plot(X(1,id),X(2,id),'ro','MarkerFace','r','MarkerEdge','k','Parent',h);
axis equal;
end



