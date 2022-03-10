function [V] = getTangentsCircle(X,M,F)
%GETTANGENTSCIRCLE for each point return its tangent model
d = 2; % dimension of the tangent model
n = size(X,2);
V = nan(d,n);
%figure; hold all;
for i = 1:n
    if(F(i)>0)
        circle = M(:,F(i));
        V(:,i) = tangentOnCircle(circle,X(:,i));
        %plot(X(1,i),X(2,i));
        %quiver(X(1,i),X(2,i),V(1,i),V(2,i));
    end
end



end

