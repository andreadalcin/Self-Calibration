function [Pcleansed, isMeaningful] = cleansePrefMatPost(R, P, delta ,kappa,b, epsiNfa)
%CLEANSEPREFMAT statistical validation step for preprocessing
% reference: Tepper and Sapiro Fast L1-NMF for multiple Parametric model
% estimation. Section 2.5
% R resiudal matrix
% P preference matrix
% kappa = locality parameter
% b = cardinality of the minimal sample set
% epsi = threhsold on false alarm
if(isempty(R)|| isempty(P))
    Pcleansed = [];
    isMeaningful = [];
    return;
end

warning('off','MATLAB:nchoosek:LargeCoefficient' );
if(nargin < 6)
    epsiNfa = 1;
end
invKappa = 1/kappa;
m = size(R,2); %number of models
isMeaningful = false(1,m);
ntest = nchoosek(m,b);

%% order clusters per cardinality

[~,clustId] = sort(sum(R<delta,1),'descend');

%%
for j = clustId
    if(sum(R(:,j)<=delta)<=b)
        isMeaningful(j) = false;
    else
        nFalseAlarm = nfa(R,j,delta,kappa, invKappa,b, ntest);
        isMeaningful(j) = (nFalseAlarm< epsiNfa);
        if(isMeaningful(j))
            R(R(:,j)<=delta,:) = Inf;
        end
    end
end

Pcleansed = P(:,isMeaningful);
warning('on','MATLAB:nchoosek:LargeCoefficient' );
end


function y = nfa(R, j, delta, kappa, invKappa,b, ntest)
% R matrix of residuals
% j index of a model in the preference matrix
% delta is the inlier threhsold
% invKappa is the inverse of the locality parameter kappa
% b cardinality of the minimal sample set
% ntest = nchoosek(m,b) number of points and minimal sample set
Cs_d = R(:,j)< delta;
m_d = sum(Cs_d);
Cs_kd = R(:,j) < kappa*delta;
m_kd = sum(Cs_kd);
y = ntest* binomialTail(m_kd -b, m_d-b, invKappa);
end




function b =  binomialTail(l,k,p)
if(p*l< k && k < l)
    r = k/l;
    b = exp(-l*(r*log(r/p)+ (1-r)*log((1-r)/(1-p))));
    return;
else
    b = 0;
    for i = k:l
        b = b+nchoosek(l,i)*(p^i)*(1-p)^(l-i);
    end
end

% approximate computation
% from https://link.springer.com/content/pdf/10.1007/978-0-387-74378-3_4.pdf

end