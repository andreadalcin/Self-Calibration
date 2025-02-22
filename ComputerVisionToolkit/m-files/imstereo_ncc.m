function [dmap,cost]=imstereo_ncc(imL,imR,drange,ws)
    %IMSTEREO_NCC Stereo block-matching with NCC
    
    dmin=drange(1); dmax=drange(2);
    s1 = filter2(ones(ws),imL.^2);
    ncc = ones([size(imL),dmax-dmin+1])*Inf;
    for d=0:dmax-dmin
        imR_d = circshift(imR,[0, -(dmin+d)]);
        imR_d(:, end-(dmin+d):end) = 0;
        prod =  imL.*imR_d; % product
        s2 = filter2(ones(ws),imR_d.^2);
        ncc(:,:,d+1) = 1-(filter2(ones(ws),prod)./sqrt(s1.*s2));
    end
    [cost,dmap]=min(ncc,[],3);
    dmap=dmap+dmin-1;
end

    
    
    
