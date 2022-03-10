function par = lineFromTangent(v,x)
%LINEFROMTANGENT Given a tangent vector and a point compute the equation of
%the line. par are [a;b;c] in the ax +by +c = 0 equation.

tol = 1e-8;
par = nan(3,1);
if( abs(v(1))<tol)
    %vertical line
    par(1) = 1;
    par(2) = 0;
    par(3) = x(1);
    return;
end

m = v(2)/v(1);
q = x(2)-m*x(1);

par(1) = m;
par(2) = -1;
par(3) = q;
par = par./norm(par(1:2));


end

