function Y = sampler_homof( X, opts)
%SAMPLER function for sampling hypotheses for MCT - T-Linkage.
% 
% opts.sampling indicates the policy of sampling:
%
% 'uniform' sampling, mss are randomly extracted.
%
% 'localized' sampling, euclidean distance is used to promote the
% extraction of closest mss. 'opts.quantile'is a parameter the controls
% probability of extracting close points. The smaller, the less the
% probability of extracting points more distant than quantile.
%
% 'nearest' as localized sampling promote the extraction of close mss.
% 'min_neigh' specifies the number of neighborood in which mss are
% extracted.
%
% opts.robust = 'x84' a robust refitting is performed once a new model
% hypotesis is instantiated. 'off' to turn off this feature.
%
% opts.geo = 1 if geometric distances are used to compute resiudals.
% 
% opts.voting specifies the different voting functions used to express
% preferences of points.
%
%  opts.voting = 'binary' opts.epsi is used as cutoff to express binary
%  voting as in J-Linkage.
%
% opts.voting = 'exp' exponential voting
% 
% opts.voting = 'gauss' gaussian voting, as exp but with 0 derivative in 0
% and a more gentle cutoff
% 
% opts.voting = 'guassx84' as 'gauss' but the influence of the threshold is
% reduced as a robust refinement is used.
%
% Please cite 
% L. Magri, A. Fusiello; Fitting Multiple Heterogeneous Models by Multi-class Cascaded T-linkage CVPR, 2019. 

if(~isfield(opts,'sampling'))
    opts.sampling = 'uniform';
    fprintf('\t uniform sampling is adopted\n');
end
if(~isfield(opts,'robust') )
    opts.robust = 'off';
    fprintf('\t robust refitting is turned off\n')
end
if(~isfield(opts,'geo'))
    opts.geo = 0;
end
if(~isfield(opts,'voting'))
    opts.voting = 'binary';
    fprintf('\t binary voting is adopted\n');
end

if(opts.geo==1)
    %fprintf('\t geometric distances are adopted\n');
end


epsi      = opts.epsi;
model     = opts.model;
sampling  = opts.sampling;
m         = opts.m;
robust    = opts.robust;
voting    = opts.voting;

n = size(X,2);

switch sampling
    case 'uniform'
        sampling_id = 0;
    case 'localized'
        if( ~isfield(opts,'distance'))
            D = squareform(pdist(X','euclidean'));
        else
            D = opts.distance;
        end
        if( ~isfield(opts,'quantile'))
            qnt = 0.65;
        else
            qnt = opts.quantile;
        end
        sigma = quantile(D(:),qnt)/4;
        sampling_id = 1;
    case 'nearest'
        tree = KDTreeSearcher(X');
        
        num_neigh = min( opts.num_neigh, size(X,2));
        
        sampling_id = 2;
    otherwise
        sampling = 'uniform';
        sampling_id = 0;
        warning('unifrom sampling is adopted');
end


switch robust
    case 'off'
        robust_id = 0;
    case 'x84'
        robust_id = 1;
    otherwise
        robust_id = 0;
end



switch voting
    case 'binary'
        voting_id = 1;
    case 'exp'
        voting_id = 2;
    case 'gauss'
        voting_id = 3;
    case 'gaussx84'
        voting_id = 4;
    otherwise
        voting_id = 0;
end

R = nan(n,m);


switch opts.model
    
    
    case 'line'
        cardmss = 2;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
                % create model
                h1 = fit_line(mss);
                flg_degen = is_invalid_line(h1);
                cont_degen = cont_degen+1;
                
            end
            
            if(flg_degen == 1)
                inds = nan(3,1);
                h1 = rand(3,1);
            end
            
            S(:,j) = inds;
            
            
            % compute residuals
            d1 = abs(h1(1).*X(1,:) + h1(2)*X(2,:) + h1(3));
            
            % refit model
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            
            card_cs1 = size(inliers,2);
            if(card_cs1==0)
                keyboard
            end
            
            
            if(card_cs1 > cardmss)
                h2 = fit_line(inliers);
                d2 = abs(h2(1).*X(1,:) + h2(2)*X(2,:) + h2(3));
                card_cs2 = sum(d2<=epsi);
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
        
        
        
    case 'line_parallel'
        L = opts.line;
        % compute angular coefficient
        m_line = -L(1)/L(2);
        
        cardmss = 1;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            %if(strcmp(sampling,'uniform'))
            if( sampling_id == 0)
                inds = randsample(n, cardmss,false);
                %elseif(strcmp(sampling,'nearest'))
            elseif(sampling_id == 2)
                seed = randsample(n, 1);
                bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                inds = datasample(bucket, cardmss,'Replace',false);
                %elseif(strcmp(sampling,'localized'))
            elseif(sampling_id == 1)
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
            h1 = fit_line_parallel( mss, m_line );
            
            % compute residuals
            %d1 = abs(h1(1).*X(1,:) + h1(2)*X(2,:) + h1(3));
            d1 = res_line(X,h1);
            %             % refit model
            %             %if(strcmp(robust,'x84'))
            %             if(robust_id)
            %                 epsi_x84 = x84(d1(d1<=epsi),cardmss);
            %                 inliers = X(:,d1<=min(epsi,epsi_x84));
            %             else
            %                 inliers = X(:,d1<=epsi);
            %             end
            %
            %             card_cs1 = size(inliers,2);
            %             assert(card_cs1~=0)
            %
            %
            %             if(card_cs1 > cardmss)
            %                 h2 = fit_line(inliers);
            %                 d2 = abs(h2(1).*X(1,:) + h2(2)*X(2,:) + h2(3));
            %                 card_cs2 = sum(d2<=epsi);
            %             else
            %                 card_cs2 = 0;
            %             end
            %             % store sampling result
            %
            %
            %             if(card_cs2 >= card_cs1)
            %                 H(:,j) = h2;
            %                 R(:,j) = d2;
            %             else
            %                 H(:,j) = h1;
            %                 R(:,j) = d1;
            %             end
            %         end
            H(:,j) = h1;
            R(:,j) = d1;
        end
        P = double(R<=epsi);
        
    case 'line_from_circle'
        center = opts.circle(1:2);
        % compute angular coefficient
        
        
        cardmss = 1;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            %if(strcmp(sampling,'uniform'))
            if( sampling_id == 0)
                inds = randsample(n, cardmss,false);
                %elseif(strcmp(sampling,'nearest'))
            elseif(sampling_id == 2)
                seed = randsample(n, 1);
                bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                inds = datasample(bucket, cardmss,'Replace',false);
                %elseif(strcmp(sampling,'localized'))
            elseif(sampling_id == 1)
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
            h1 = fit_line_from_circle( mss, center);
            
            % compute residuals
            %d1 = abs(h1(1).*X(1,:) + h1(2)*X(2,:) + h1(3));
            d1 = res_line(X,h1);
            %             % refit model
            %             %if(strcmp(robust,'x84'))
            %             if(robust_id)
            %                 epsi_x84 = x84(d1(d1<=epsi),cardmss);
            %                 inliers = X(:,d1<=min(epsi,epsi_x84));
            %             else
            %                 inliers = X(:,d1<=epsi);
            %             end
            %
            %             card_cs1 = size(inliers,2);
            %             assert(card_cs1~=0)
            %
            %
            %             if(card_cs1 > cardmss)
            %                 h2 = fit_line(inliers);
            %                 d2 = abs(h2(1).*X(1,:) + h2(2)*X(2,:) + h2(3));
            %                 card_cs2 = sum(d2<=epsi);
            %             else
            %                 card_cs2 = 0;
            %             end
            %             % store sampling result
            %
            %
            %             if(card_cs2 >= card_cs1)
            %                 H(:,j) = h2;
            %                 R(:,j) = d2;
            %             else
            %                 H(:,j) = h1;
            %                 R(:,j) = d1;
            %             end
            %         end
            H(:,j) = h1;
            R(:,j) = d1;
        end
        P = double(R<=epsi);
        
    case 'circle'
        cardmss = 3;
        H = nan(3,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if(sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
            if(flg_degen == 1)
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_circle_taubin(mss);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = abs( sqrt(sum((X-repmat(h1(1:2),1,n)).^2,1))-h1(3));
            
            % refit model
            
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                h2 = fit_circle_lm(inliers,h1);
                d2 = abs( sqrt(sum((X-repmat(h2(1:2),1,n)).^2,1))-h2(3));
                card_cs2 = sum(d2<=epsi);
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
        
    case 'parabola'
        cardmss = 3;
        S = nan(cardmss,m);
        H = nan(3,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
                % create model
               
                flg_degen = iscolinear(mss(:,1),mss(:,2),mss(:,3));
                cont_degen = cont_degen+1;
                %h1 = fit_parabola(mss);
              
                
            end
            
            if(flg_degen == 1)
                disp('nooo')
                inds = nan(3,1);
                h1 = rand(3,1);
            else
                h1 = fit_parabola(mss);
            end
            
            S(:,j) = inds;
            
            
            % compute residuals
            d1 = res_parabola(X,h1);
%             figure; hold all;
%             scatter(X(1,:),X(2,:),d1)
%             plot(mss(1,:),mss(2,:),'r*')
%             x1 = linspace(min(X(1,:)),max(X(1,:)));
%             y1 = polyval(h1,x1);
%             plot(x1,y1);
            % refit model
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            
            card_cs1 = size(inliers,2);
            if(card_cs1==0)
                %keyboard
            end
            
            
            if(card_cs1 > cardmss)
                h2 = fit_parabola(inliers);
                d2 = res_parabola(X,h2);
                card_cs2 = sum(d2<=epsi);
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
        
        
    case 'fundamental'
        %disp('fundamental matrix fitting')
        
        %
        cardmss = 8;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
                f_tentative = fit_fm(mss);
                flg_degen = is_fundamental_degen( mss, f_tentative , cardmss );
                %                 if(flg_degen==1)
                %                     disp('degenerate')
                %                 end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = f_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo==1)
                d1 = res_fm_geo(X, h1);
            else
                d1 = res_fm(X, h1);
            end
            % refit model
            
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                
                %h2 = fit_fm_torr(inliers,h1);
                h2 = fit_fm(inliers);
                
                if(opts.geo==1)
                    d2 = res_fm_geo(X, h2);
                else
                    d2 = res_fm(X, h2);
                end
                
                card_cs2 = sum(d2<=epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
            if((card_cs2 >= card_cs1) && card_cs2>0)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
            
        end
        
        
        
        %% %%%%
    case 'affine_fundamental'
        disp('affine fundamental matrix fitting')
        cardmss = 4;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
                f_tentative = fit_fm_affine(mss);
                flg_degen = 0; %is_homography_degen( mss, f_tentative , cardmss );
                %                 if(flg_degen==1)
                %                     disp('degenerate')
                %                 end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = f_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo==1)
                d1 = res_fm_geo(X, h1);
            else
                d1 = res_fm(X, h1);
            end
            % refit model
            
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                
                %h2 = fit_fm_torr(inliers,h1);
                
                h2 = fit_fm(inliers);
                
                if(opts.geo==1)
                    d2 = res_fm_geo(X, h2);
                else
                    d2 = res_fm(X, h2);
                end
                card_cs2 = sum(d2<=epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
            if((card_cs2 >= card_cs1) && card_cs2>0)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
            
        end
        
        
        
        %% %%%%
    case 'homography'
        cardmss = 4;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 50;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
                if(~validateMSS_homography(X, inds))
                    flg_degen = 1 ;
                else
                    h_tentative = fit_homography(mss);
                    flg_degen = is_homography_degen(X, h_tentative, inds);
                end
                %if(flg_degen==1)
                %    disp('degenerate')
                %end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                disp('rand')
                %inds = randsample(n, cardmss,false);
                % mss = X(:,inds);
                % h1 = fit_homography(mss);
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                h1 = h_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            if(opts.geo == 1)
                d1 = res_homography_geo(X, h1);
            else
                d1 = res_homography(X,h1);
            end
            % refit model
            
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                if(card_cs1>3*cardmss)
                    %h2 = fit_homography_nonlin(inliers,reshape(h1,[3,3]));
                    h2 = fit_homography(inliers);
                    h2 = h2(:);
                    
                else
                    h2 = fit_homography(inliers);
                end
                
                if(rcond(reshape(h2,[3,3]))<1e-15)
                    h2 = h1;
                end
                if(opts.geo==1)
                    d2 = res_homography_geo(X, h2);
                else
                    d2 = res_homography(X, h2);
                end
                card_cs2 = sum(d2<=epsi);
                
            else
                card_cs2 = 0;
            end
            % store sampling result
            if((card_cs2 >= card_cs1) && card_cs2>0)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
            
        end
        
        
        
        %%
        
        
    case 'homography_from_fund'
        
        f = opts.fund;
        F = reshape(f,[3,3]);
        e2 = epipole(F');%null(F');
        A = star(e2)*F;
        
        cardmss = 3;
        H = nan(9,m);
        S = nan(cardmss,m);
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 20;
            while(flg_degen && cont_degen<=max_cont_degen)
                
                %instantiate mss
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    % uniform sampling
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                    %elseif(strcmp(sampling,'localized'))
                elseif(sampling_id == 1)
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
                x1 = X(1:3,inds);
                x2 = X(4:6,inds);
                flg1 = iscolinear(x1(:,1),x1(:,2),x1(:,3),'h')|| iscolinear(x2(:,1),x2(:,2),x2(:,3),'h');
                if(~flg1)
                    
                    h_tentative = fit_homography_from_fund(mss,A,e2);
                    
                    flg_degen = is_homography_from_fund_degen(X, h_tentative, F, inds);
                end
                
                %if(flg_degen==1)
                %    disp('degenerate')
                %end
                
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(cardmss,1);
                h1 = rand(9,1);
            else
                
                h1 = h_tentative;
            end
            S(:,j) = inds;
            
            % compute residuals
            
            d1 = res_homography(X, h1);
            
            % refit model
            
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            card_cs1 = size(inliers,2);
            
            
            if(card_cs1 > cardmss)
                if(card_cs1>4*cardmss)
                    %h2 = fit_homography_nonlin(inliers,reshape(h1,[3,3]));
                    h2 = fit_homography(inliers);
                    h2 = h2(:);
                else
                    h2 = fit_homography(inliers);
                end
                if(rcond(reshape(h2,[3,3]))<1e-15)
                    h2 = h1;
                end
                d2 = res_homography(X, h2);
                card_cs2 = sum(d2<=epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
            if((card_cs2 >= card_cs1) && card_cs2>0)
                H(:,j) = h2;
                R(:,j) = d2;
            else
                H(:,j) = h1;
                R(:,j) = d1;
            end
            
        end
    case 'subspace'
        
        delta = opts.dim_subspace;
        cardmss = delta;
        f =size(X,1);
        S = nan(cardmss,m);
        H = nan(f*f,m);
        % main loopp
        for j = 1:m
            % instantiate mss
            flg_degen = true;
            cont_degen = 0;
            max_cont_degen = 10;
            while(flg_degen && cont_degen<=max_cont_degen)
                %if(strcmp(sampling,'uniform'))
                if( sampling_id == 0)
                    inds = randsample(n, cardmss,false);
                    %elseif(strcmp(sampling,'nearest'))
                elseif(sampling_id == 2)
                    seed = randsample(n, 1);
                    bucket = knnsearch(tree,X(:,seed)','k',num_neigh);
                    inds = datasample(bucket, cardmss,'Replace',false);
                elseif(sampling_id == 1)
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
                % create model
                h1 = fit_subspace(mss,delta);
                flg_degen = is_invalid_subspace(mss,delta);
                cont_degen = cont_degen+1;
                
            end
            
            if(flg_degen == 1)
                disp('nooo')
                inds = nan(3,1);
                h1 = rand(3,1);
            end
            
            S(:,j) = inds;
            
            
            % compute residuals
            d1 = res_subspace(X,h1);
            
            % refit model
            %if(strcmp(robust,'x84'))
            if(robust_id)
                epsi_x84 = x84(d1(d1<=epsi),cardmss);
                inliers = X(:,d1<=min(epsi,epsi_x84));
            else
                inliers = X(:,d1<=epsi);
            end
            
            card_cs1 = size(inliers,2);
            if(card_cs1==0)
                keyboard
            end
            
            
            if(card_cs1 > cardmss)
                h2 = fit_subspace(inliers,delta);
                d2 = res_subspace(X,h2);
                card_cs2 = sum(d2<=epsi);
            else
                card_cs2 = 0;
            end
            % store sampling result
            
            
            if(card_cs2 >= card_cs1)
                H(:,j) = h2.P(:);
                R(:,j) = d2;
            else
                H(:,j) = h1.P(:);
                R(:,j) = d1;
            end
        end
        
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
        
    case 'cylinder'
        disp('refit not supported')
        H = nan(7,m);
        cardmss = 2;
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
                flg_degen = is_cylinder_degen(mss);
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(2,1);
                h1 = rand(7,1);
            else
                h1 = fit_pc_cylinder(mss);
                w0 = h1(4:6); % direction of the cylinder;
                h1 = convertToFiniteCylinder(h1,X);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_pc_cylinder(X,h1);
            
            % refit model
            
            if(strcmp(robust,'x84'))
                epsi_x84 = x84(d1(d1<epsi),cardmss);
                inliers = X(:,d1<min(epsi,epsi_x84));
            else
                inliers = X(:,d1<epsi);
            end
            card_cs1 = size(inliers,2);
            if(card_cs1<1)
                %keyboard
            end
            
            if(card_cs1 > 20) % refit not supported
                h2 = fit_pc_cylinder_ls(inliers,w0);
                h2 = convertToFiniteCylinder(h2,X);
                d2 = res_pc_cylinder(X,h2);
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
    case 'cylinder_axis'
        disp('cylinder with fixed direction')
        w0 = opts.axis;
        H = nan(7,m);
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
                flg_degen = is_cylinder_degen(mss);
                cont_degen = cont_degen+1;
            end
            
            % create model
            if(flg_degen == 1)
                inds = nan(2,1);
                h1 = rand(7,1);
            else
                h1 = fit_pc_cylinder_ls_axis(mss,w0);
                h1 = convertToFiniteCylinder(h1,X);
            end
            S(:,j) = inds;
            
            % compute residuals
            d1 = res_pc_cylinder(X,h1);
            
            % refit model
            
            if(strcmp(robust,'x84'))
                epsi_x84 = x84(d1(d1<epsi),cardmss);
                inliers = X(:,d1<min(epsi,epsi_x84));
            else
                inliers = X(:,d1<epsi);
            end
            card_cs1 = size(inliers,2);
            if(card_cs1<1)
                %keyboard
            end
            
            if(card_cs1 > 20) % refit not supported
                h2 = fit_pc_cylinder_ls_axis(inliers,w0);
                h2 = convertToFiniteCylinder(h2,X);
                d2 = res_pc_cylinder(X,h2);
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
        
        
    otherwise
        error('the model is not defined: sampler supports ''line'' and ''circle'', ''homography'' and ''fundamental''.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch voting_id
    case 1
        P = double(R<=epsi);
    case 2
        % exponential voting
        P = zeros(size(R));
        tau=epsi/5;
        I = R<epsi;
        P(I) = exp(-R(I)./tau);
    case 3
        P = zeros(size(R));
        heta = 0.05; % minumum preference for residual equal to epsilon
        sigma2 = -epsi^2/log(heta);
        I = R<epsi;
        P(I) = exp(-(R(I).^2)./(sigma2));
    case 4
        % compute threshold per column
        
        theta =  2.5;
        location = nanmedian(R,1);
        scale = theta/0.6745 * nanmedian(abs(R-location),1);
        thresh = min(scale + location, epsi);
        
        % vote
        heta = 0.05; % minumum preference for residual equal to epsilon
        sigma2 =  -thresh.^2/log(heta);
        V = exp(-(R.^2)./(sigma2));
        
        P = zeros(size(R));
        I = R< thresh;
        P(I) = V(I);
        %%
end




Y.S = S;
Y.H = H;
Y.R = R;
Y.P = P;





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

