function [ok, ms, msOutput] = isMergeableGricFlatsMixed(X, L, i, j, lambda1 , lambda2, sigma,dimFlats)
% Check if two clusters A and B can be merged.
o_verbose = false;
[okl, msf1, msOutf1] = isMergeableGricFlat(X, L, i, j, lambda1 , lambda2, sigma, dimFlats(1));
[okc, msf2, msOutf2] = isMergeableGricFlat(X, L, i, j, lambda1 , lambda2, sigma, dimFlats(2));
[okc, msf3, msOutf3] = isMergeableGricFlat(X, L, i, j, lambda1 , lambda2, sigma, dimFlats(3));
%%
scores = [msf1.gric.before, msf1.gric.after,...
    msf2.gric.before, msf2.gric.after,...
    msf3.gric.before, msf3.gric.after];
[~, ind] = min(scores);

if(o_verbose)
    f = [msf1.fidelity.before, msf1.fidelity.after,...
        msf2.fidelity.before, msf2.fidelity.after...
        msf3.fidelity.before, msf3.fidelity.after];
    c = [msf1.complexity.before, msf1.complexity.after, msf2.complexity.before, msf2.complexity.after]./lambda2;
    disp('--------------------------------------')
    disp([scores;f;c]);
    disp('--------------------------------------')
end
%%  keep the result with the minimum model selection score
switch(ind)
    case 1
        % do not merge clusters: seprate flats of lower dimension win
        ok = false;
        ms = msf1;
        msOutput = msOutf1;
    case 2
        % do merge clusters with flat
        ok = true;
        ms = msf1;
        msOutput = msOutf1;
        assert(okl);
        
    case 3
        % do not merge clustes: separate flats of bigger dimension win
        ok = false;
        ms = msf2;
        msOutput = msOutf2;
    case 4
        % do merge clusters with flat
        ok = true;
        ms = msf2;
        msOutput = msOutf2;
        assert(okc);
    case 5
        % do not merge clustes: separate flats of bigger dimension win
        ok = false;
        ms = msf3;
        msOutput = msOutf3;
    case 6
        % do merge clusters with flat
        ok = true;
        ms = msf3;
        msOutput = msOutf3;
        assert(okc);
end




end
