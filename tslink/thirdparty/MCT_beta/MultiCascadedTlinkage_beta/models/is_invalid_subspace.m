function [b] = is_invalid_subspace(mss,delta)
%IS_INVALID_SUBSPACE 
b = rank(mss)<delta;

end

