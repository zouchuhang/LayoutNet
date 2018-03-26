function [pMean, pVal, estBias, nv] = evaluateProbabilityEstimate(pEst, pTrue, step)

pVal = [step/2:step:(1-step/2)];

nsteps = numel(pVal);

estBias = 0;

for v = 1:numel(pVal)
    v1 = pVal(v)-step/2;
    v2 = min(v1+step, 1);
    
    ind = find((pEst >= v1) & (pEst <= v2));
    nv(v) = numel(ind);
    if numel(ind)>0
        pMean(v) = mean(pTrue(ind));
        estBias = estBias + numel(ind)/numel(pEst)*(pMean(v)-pVal(v));
    else
        pMean(v) = 0;
    end
end


