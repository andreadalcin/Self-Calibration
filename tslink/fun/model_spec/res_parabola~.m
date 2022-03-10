function resi = res_parabola(X,h)
m = size(h,2);
n = size(X,2);
resi = nan(n,m);
for j = 1:m
    H = h(:,j);
    for ii =1:n
        P = X(:,ii);
        
        gamma = H(3)-P(2);
        % derivative
        f = [4*H(1)^2, 6*H(1)*H(2), 2*( 1 +H(2)^2 + 2*H(1)*gamma), -2*P(1) + 2*H(2)*gamma ];
        u = roots(f);
        
        
        k = numel(u); % number of solutions
        Q = nan(2,k); % possible extrema points on the parabola
        d = Inf(1,k);
        
        for i=1:k
            if(isreal(u(i)))
                Q(:,i) = [u(i);H(1)*u(i)^2+H(2)*u(i)+H(3)];
                d   = norm(Q(:,i)-P);
            end
        end
        [mindist,~] = min(d);
        resi(ii,j) = mindist;
    end
end
end