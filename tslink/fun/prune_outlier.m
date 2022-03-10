function [F] = prune_outlier(X,C,modelType,epsi ,kappa, epsiNfa)
%PRUNE_OUTLIER Summary of this function goes here
%   Detailed explanation goes here
F = C;

o_debugPlot = false;

if(strcmp(modelType,'lc'))
    cardmss_l= 2;
    cardmss_c = 3;
    % line and circle
    Ml = recover_line(X,F);
    Rl = res_line(X,Ml);
    Pl = resiToP(Rl,epsi);
    [~, isMeaningful_l] = cleansePrefMatPost(Rl, Pl, epsi ,kappa,cardmss_l, epsiNfa);
    if(o_debugPlot)
        figure; hold all;
        subplot(1,3,1);
        hold all;
        for i = 1:max(C)
            if(isMeaningful_l(i))
                scatter(X(1,C==i),X(2,C==i));
            end
        end
    end
    Mc = recover_circle(X,F);
    Rc = res_circle(X,Mc);
    Pc = resiToP(Rc,epsi);
    [~, isMeaningful_c] = cleansePrefMatPost(Rc, Pc, epsi ,kappa,cardmss_c, epsiNfa);
    
    if(o_debugPlot)
        subplot(1,3,2);
         hold all;
        for i = 1:max(C)
            if(isMeaningful_c(i))
                scatter(X(1,C==i),X(2,C==i));
            end
        end
    end
    
    
    
    isMeaningful = isMeaningful_c | isMeaningful_l;
    
    if(o_debugPlot)
        subplot(1,3,3);
         hold all;
        for i = 1:max(C)
            if(isMeaningful(i))
                scatter(X(1,C==i),X(2,C==i));
            end
        end
    end
    
elseif(strcmp(modelType,'lcp'))
    cardmss_l= 2;
    cardmss_c = 3;
    cardmss_p = 3;
    % line circle and parabolas
    Ml = recover_line(X,F);
    Rl = res_line(X,Ml);
    Pl = resiToP(Rl,epsi);
    [~, isMeaningful_l] = cleansePrefMatPost(Rl, Pl, epsi ,kappa,cardmss_l, epsiNfa);
    
    Mc = recover_circle(X,F);
    Rc = res_circle(X,Mc);
    Pc = resiToP(Rc,epsi);
    [~, isMeaningful_c] = cleansePrefMatPost(Rc, Pc, epsi ,kappa,cardmss_c, epsiNfa);
    
    Mp = recover_parabola(X,F);
    Rp = res_parabola(X,Mp);
    Pp = resiToP(Rc,epsi);
    [~, isMeaningful_p] = cleansePrefMatPost(Rp, Pp, epsi ,kappa,cardmss_p, epsiNfa);
    
    isMeaningful = isMeaningful_c | isMeaningful_l | isMeaningful_p;
    
end


for j = 1:max(F)
    if(~isMeaningful(j))
        F(F==j) = 0;
    end
end
F(F~=0)=grp2idx(F(F~=0));


end

