function guiPrefExplorer(tree,X,S,radius,useTani, picked, J,T,cmap)
f = figure;
title('Preference Explorer')
% initialization
useTani = false;
T = [];
J = [];
picked = false(size(X,2),1);
cmap = brewermap(300,'YlGnBu');
radiusNeigh = 1;
s = cell(1,numel(picked)); % handle to points and preferred models
if(strcmp(S.model,'circle') )
    useCircle = true;
else
    useCircle = false;
end
%% ----------------------------------------------------------------------
buttonShow = uicontrol;
buttonShow.Position = [20 20 60 20];
buttonShow.String = 'Show';
buttonShow.Callback = @buttonShowPushed;

    function  buttonShowPushed(src,event)
        
        scatter(X(1,:),X(2,:),'k','filled');
        axis equal;
        
    end
%% ----------------------------------------------------------------------
buttonPick = uicontrol;
buttonPick.Position = [90 20 60 20];
buttonPick.String = 'Pick';
buttonPick.Callback = @buttonPickPushed;

    function  buttonPickPushed(src,event)
        fprintf('Right click to stop picking points\n');
        [picked,s] = pickPointAndPrefModels(tree,X,S,radius,picked,s,useCircle);
    end

%% ----------------------------------------------------------------------
buttonIntersect = uicontrol;
buttonIntersect.Position = [160 20 60 20];
buttonIntersect.String = 'Intersect';
buttonIntersect.Callback = @buttonIntersectPushed;

    function buttonIntersectPushed(src,event)
        if(sum(picked)>1)
            %cla
             displayIntersecPreferredModels(picked, X, S, radius,cmap, useCircle);
        else
            disp('Not enough points selected');
        end
    end

%% ----------------------------------------------------------------------
buttonDistance = uicontrol;
buttonDistance.Position = [230 20 60 20];
buttonDistance.String = 'Distance';
buttonDistance.Callback = @buttonTanimotoPushed;

    function buttonTanimotoPushed(src,event)
        if(useTani)
            disp('Use Tanimoto distances');
            if(isempty(T))
                n = size(X,2);
                s0 = S.P*S.P'; % inner product
                n0 = diag(s0); %norm
                m0 = repmat(n0,1,n);
                d0 = m0 + m0' - s0;
                T =  1 - s0./d0;
            end
            
        else
            disp('Use Jaccard distances');
            if(isempty(J))
                J =  squareform(pdist(S.P,'jaccard'));
            end
        end
    end
%% ----------------------------------------------------------------------
buttonNiegh = uicontrol;
buttonNiegh.Position = [300 20 60 20];
buttonNiegh.String = 'Neigh';
buttonNiegh.Callback = @buttonNeighPushed;

    function buttonNeighPushed(src,event)
        if(sum(picked)==0)
            disp('Pick point before');
            return
        end
        
        if(isempty(J))
            disp('Calculate distance before.');
            return;
        end
        if(sum(picked)>2)
            disp('Be patient :)')
        end
        fprintf('Computing and rendering distance map...')
        for i = 1:size(X,2)
            if(picked(i))
                dists = J(i,:);
                displayPrefNeigh(X,dists,radiusNeigh,1.2*radius,cmap);
            end
        end
        fprintf('done\n');
    end




%% ----------------------------------------------------------------------
sliderRadius = uicontrol('Style','Slider','min',0,'max',1,'value',1);
sliderRadius.Position = [370 20 60 20];
sliderRadius.String = 'Radius';
sliderRadius.Callback = @sliderRadiusPushed;

    function sliderRadiusPushed(src,event)
        num = get(sliderRadius,'value');
        fprintf('Value of the radius is %.2d\n',num);
        radiusNeigh = num;
    end
%% ----------------------------------------------------------------------
buttonReset = uicontrol;
buttonReset.Position = [440 20 60 20];
buttonReset.String = 'Reset';
buttonReset.Callback = @buttonResetPushed;

    function buttonResetPushed(src,event)
        picked = false(size(X,2),1);
        cla;
    end
end