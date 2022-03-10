function Y = sampler_motion( X, opts)
%SAMPLER
% 2do dbscan
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
    case 'uniform'
        sampling_id = 0;
    case 'localized'
        D = squareform(pdist(X','euclidean'));
        sigma = max(D(:))/4;
        sampling_id = 1;
    case 'nearest'
        tree = KDTreeSearcher(X');
        
        num_neigh = opts.num_neigh;
        
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
                disp('nooo')
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
        
        P = double(R<=epsi);
        
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
        
        P = double(R<=epsi);
        
    case 'fundamental'
        disp('fundamental matrix fitting')
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
            
            d1 = res_fm(X, h1);
            
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
                
                
                d2 = res_fm(X, h2);
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
        
        P = double(R<=epsi);
        
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
            
            d1 = res_fm(X, h1);
            
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
                
                
                d2 = res_fm(X, h2);
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
        
        P = double(R<=epsi);
        
        %% %%%%
    case 'homography'
        cardmss = 4;
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
                h_tentative = fit_homography(mss);
                flg_degen = is_homography_degen(X, h_tentative, inds);
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
        
        P = double(R<epsi);
        
        %%
        
        
    case 'homography_from_fund'
        
        f = opts.fund;
        F = reshape(f,[3,3]);
        e2 = null(F');
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
        
        P = double(R<=epsi);
        
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

