function [] = displayIntersecPreferredModels(active, X, S, radius,cmap, useCircle)
%DISPLAYINERSECTIONPREFERREDMODELS
h = hggroup('Parent', gca);
prefIntersect = prod(S.P(active,:),1);
% set colors, if available use colormap
if(nargin<5)
    if(sum(prefIntersect)==0)
        colr = 'k';
    else
        colr = 'c';
    end
else
    if(sum(prefIntersect)==0)
        colr = cmap(floor(size(cmap,1)/2),:);
    else
        colr = [14,71,158]./255;cmap(floor(size(cmap,1)/2),:);cmap(end,:);
    end
end
% plot preferred models

if(useCircle)
    for id = 1:size(X,2)
        if(active(id))
            point = X(:,id);
            drawCircle(point',radius,h);
            for j = 1:numel(prefIntersect)
                if(prefIntersect(j)>0)
                    
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
                        drawArcCircle(S.H(1:2,j),S.H(3,j), thetaStart, thetaEnd,h, colr, 2)
                    end
                end
            end
        end
    end
    
else
    
    for id = 1:size(X,2)
        if(active(id))
            point = X(:,id);
            drawCircle(point',radius,h);
            for j = 1:numel(prefIntersect)
                if(prefIntersect(j)>0)
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
                    
                    [xout,yout] = linecirc(slope,intercpt,point(1),point(2),radius);
                    line(xout,yout,'Color',colr,'LineWidth',2,'Parent',h)
                end
            end
            plot(X(1,id),X(2,id),'bo','MarkerFace',colr,'MarkerEdge','w','Parent',h);
        end
    end
end
end

