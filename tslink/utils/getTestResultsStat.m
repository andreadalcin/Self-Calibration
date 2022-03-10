function [testResults] = getTestResultsStat(testResults)
%GETTESTRESULTSSTAT Summary of this function goes here
%   Detailed explanation goes here
nAlgo = numel(testResults);
for i = 1:nAlgo
    me = testResults(i).me;
    maxIter = size(me,3);
    maxSeq = size(me,1);
    assert(maxSeq ==1);
    s =1;
    if(maxIter>0)
        % compute quartile ranges
        
        meQrtl1 = quantile(me(s,:,:),0.25,3);
        meQrtl3 = quantile(me(s,:,:),0.75,3);
        testResults(i).qrt.l1 = meQrtl1;
        testResults(i).qrt.l3 = meQrtl3;
        testResults(i).qrtRegion = [meQrtl1,fliplr(meQrtl3)];
        
    end
    
    testResults(i).median = median(me(s,:,:),3);
    testResults(i).mean = median(me(s,:,:),3);
    testResults(i).min = min(me(s,:,:),[],3);
    testResults(i).std = std(me(s,:,:),0,3);
    testResults(i).max = max(me(s,:,:),[],3);
    [~, bestE] = min(mean(me(s,:,:),3));
    testResults(i).bestEpsi =  testResults(i).epsisVec(bestE);
    [~,bestIter] = min(me(s,bestE,:));
    [~,worstIter] = max(me(s,bestE,:));
    testResults(i).bestMinMeanClust = testResults(i).clust{s,bestE,bestIter};
    testResults(i).worstMinMeanClust = testResults(i).clust{s,bestE,worstIter};
end
end

