function Y = sampler_plane( X, opts)
%SAMPLER
% Author: Luca Magri
if(~isfield(opts,'sampling'))
    opts.sampling = 'uniform';
    disp('uniform sampling is adopted');
end
if(~isfield(opts,'robust') )
    opts.robust = 'off';
    disp('robust refitting is turned off')
end

epsi      = opts.epsi;
model     = opts.model;
sampling  = opts.sampling;
m         = opts.m;
robust    = opts.robust;


n = size(X,2);

switch sampling
    case 'localized'
        D = squareform(pdist(X','euclidean'));
        sigma = max(D(:))/4;
    case 'uniform'
    case 'nearest'
        tree = KDTreeSearcher(X');
        num_neigh = 10;
    otherwise
        sampling = 'uniform';
        warning('unifrom sampling is adopted');
        
end


R = nan(n,m);


switch opts.model
    
   
    case 'line'
        H = nan(3,m);
        cardmss = 2;
        S = nan(cardmss,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            if(strcmp(sampling,'uniform'))
                inds = randsample(n, cardmss,false);
            elseif(strcmp(sampling,'nearest'))
                seed = randsample(n, 1);
                bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                inds = datasample(bucket, cardmss);
            elseif(strcmp(sampling,'localized'))
                inds = nan(1,cardmss);
                inds(1) = randsample(n, 1);
                w = exp(-(D(inds(1),:).^2)/sigma.^2);
                w(inds(1))=0;
                for i=2:cardmss
                    inds(i) = randsample(n,1,true,w);
                    w(inds(i)) = 0;
                end
            end
            S(:,j) = inds;
            mss = X(:,inds);
            % create model
            h1 = fit_line(mss);
            
            % compute residuals
            d1 = abs(h1(1).*X(1,:) + h1(2)*X(2,:) + h1(3));
            
            % refit model
            if(strcmp(robust,'x84'))
                epsi_x84 = x84(d1(d1<epsi),1);
                inliers = X(:,d1<min(epsi,epsi_x84));
            else
                inliers = X(:,d1<epsi);
            end
            
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                h2 = fit_line(inliers);
                d2 = abs(h2(1).*X(1,:) + h2(2)*X(2,:) + h2(3));
                card_cs2 = sum(d2<epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
       
           
            if(card_cs2 >= card_cs1)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
        end
        
        P = double(R<epsi);
        
        
        
    case 'circle'
        H = nan(3,m);
        cardmss = 3;
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                if(strcmp(sampling,'uniform'))
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                elseif(strcmp(sampling,'nearest'))
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss);
                elseif(strcmp(sampling,'localized'))
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = iscolinear(mss(:,1),mss(:,2),mss(:,3));
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1);
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_circle_taubin(mss);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = abs( sqrt(sum((X-repmat(h1(1:2),1,n)).^2,1))-h1(3));
            
            % refit model
            
            if(strcmp(robust,'x84'))
                epsi_x84 = x84(d1(d1<epsi),1);
                inliers = X(:,d1<min(epsi,epsi_x84));
            else
                inliers = X(:,d1<epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                h2 = fit_circle_lm(inliers,h1);
                d2 = abs( sqrt(sum((X-repmat(h2(1:2),1,n)).^2,1))-h2(3));
                card_cs2 = sum(d2<epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
            if(card_cs2 >= card_cs1)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
            
        end
        
        P = double(R<epsi);
        
        
        
        case 'plane'
        H = nan(4,m);
        cardmss = 3;
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                if(strcmp(sampling,'uniform'))
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                elseif(strcmp(sampling,'nearest'))
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss);
                elseif(strcmp(sampling,'localized'))
                    % localized sampling
                    inds = nan(1,cardmss);
                    inds(1) = randsample(n, 1);
                    w = exp(-(D(inds(1),:).^2)/sigma.^2);
                    w(inds(1))=0;
                    for i=2:cardmss
                        inds(i) = randsample(n,1,true,w);
                        w(inds(i)) = 0;
                    end
                end
                mss = X(:,inds);
                flg_degen = iscolinear(mss(:,1),mss(:,2),mss(:,3));
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1);
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_plane(mss);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_plane(X,h1);
            
            % refit model
            
            if(strcmp(robust,'x84'))
                epsi_x84 = x84(d1(d1<epsi),1);
                inliers = X(:,d1<min(epsi,epsi_x84));
            else
                inliers = X(:,d1<epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                h2 = fit_plane(inliers);
                d2 = res_plane(X,h2);
                card_cs2 = sum(d2<epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
            if(card_cs2 >= card_cs1)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
            
        end
        
        P = double(R<epsi);
        
        
        
          
        
    otherwise
        error('the model is not defined: sampler supports ''line'' and ''circle''.')
end







Y.S = S;
Y.H = H;
Y.R = R;
Y.P = double(P);





end

%%

function [thresh,inliers] =  x84(res, n)
% n e' il numero minimo di inliers che ritorna

theta = 2.5;
% attenzione: il theta di default di X84 e' 3.5
% cosi' e' piu' restrittivo, ma in linea con la regola di selezione degli
% inliers in LMEDS

location = nanmedian(res(:));
scale = theta/0.6745 * nanmedian(abs(res(:)-location));
inliers = abs(res-location) <= scale;
thresh = scale+location;
%fprintf('x84 scale: %f \n', scale+location);

if length(inliers) < n
    % ritorna i primi n
    [~, i] = sort(res);
    inliers = i(1:min(numel(i),n));
    inliers =  sort(inliers);
    
end
end

