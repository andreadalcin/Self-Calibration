function [] = displayMeGraph(testResults, kappaSigma)
%DISPLAYMEGRAPH
if(nargin < 2)
    kappaSigma =[];
end
nAlgo = numel(testResults);
cmap = brewermap(nAlgo,'Set2');
%cmap = cmap([2,1,3,4],:);
hold all;
p = zeros(1,nAlgo);

algoId = 1:nAlgo;
for i =algoId
    epsisVec = testResults(i).epsisVec;
    p(i) = plot(epsisVec, testResults(i).median,'linewidth',5,'Color',cmap(i,:));
end

for i = algoId
    me = testResults(i).me;
    epsisVec = testResults(i).epsisVec;
    maxIter = size(me,3);
    maxSeq = size(me,1);
    assert(maxSeq ==1);
    s =1;
    if(maxIter>0)
        % display quartile ranges
        epsiRegion = [epsisVec,fliplr(epsisVec)];
        qrtRegion = testResults(i).qrtRegion;
        fill(epsiRegion,qrtRegion,cmap(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.0);
    end
    
    
    plot(epsisVec, testResults(i).min,'o','MarkerSize',10,'linewidth',2,'Color',cmap(i,:))
    plot(epsisVec, testResults(i).max,'+','MarkerSize',10,'linewidth',2,'Color',cmap(i,:))
end
for i =algoId
    epsisVec = testResults(i).epsisVec;
    p(i) = plot(epsisVec, testResults(i).median,'linewidth',5,'Color',cmap(i,:));
end
%legend(p(1:nAlgo),{testResults.name});
ylabel('ME');
xlabel('\epsilon');
if(~isempty(kappaSigma))
    xl = {};
    xt =[];
    for i = 1:numel(epsisVec)
        if(mod(i+1,2)==0)
            xl = [xl,{[num2str(kappaSigma(i)),'\sigma']}];
            xt = [xt, epsisVec(i)];
        end
    end
    
    xticks(xt)
    xticklabels(xl)
end
set(gca,'FontSize',14);
legendNames = {};
cont = 1;
for i = algoId
    legendNames{cont} = testResults(i).name;
    cont = cont +1;
end
legend(legendNames);
hold off;
end


