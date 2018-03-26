function [estBias, estVar] = testSpSegmentationEstimator(imsegs, spdata, maps)

% UNFINISHED
[data, lab] = formatData(spdata, {imsegs(:).labels});

vconf = test_boosted_dt_mc(vclassifier, data);
vconf = 1 ./ (1+exp(-vconf));
vconf = vconf ./ repmat(sum(vconf, 2), 1, size(vconf, 2));

hconf = test_boosted_dt_mc(hclassifier, data);
hconf = 1 ./ (1+exp(-hconf));
hconf = hconf ./ repmat(sum(hconf, 2), 1, size(hconf, 2));

spEst = [vconf(:, 1) repmat(vconf(:, 2), 1, 5).*hconf2 vconf(:, 3)];
spTrue = lab;

nmaps = 0;
for f = 1:numel(maps)
    nmaps = nmaps + numel(maps{f});
end

pEst = zeros(nmaps, 1);
pTrue = zeros(nmaps, 1);
% 
% True = zeros(size(valEst));
% for i = 1:numel(lab)
%     for k = 
% 
% spvlabels = {imsegs(:).vert_labels};
% sphlabels = {imsegs(:).horz_labels};
% 
% % test vertical
% [vdata, vlab] = formatData(spfeatures, spvlabels);
% vconf = test_boosted_dt_mc(vclassifier, vdata);
% vconf = 1 ./ (1+exp(-vconf));
% vconf = vconf ./ repmat(sum(vconf, 2), 1, size(vconf, 2));
% [tmp, vmax] = max(vconf, [], 2);
% verror = mean(vmax~=vlab)
% 
% vconferr = 0;
% for k = 1:numel(vlab)
%     vconferr = vconferr + 1 - vconf(k, vlab(k));
% end
% vconferr = vconferr / numel(vlab)
% 
% % test subclass
% [hdata, hlab] = formatData(spfeatures, sphlabels);
% hconf = test_boosted_dt_mc(hclassifier, hdata);
% hconf = 1 ./ (1+exp(-hconf));
% hconf = hconf ./ repmat(sum(hconf, 2), 1, size(hconf, 2));
% [tmp, hmax] = max(hconf, [], 2);
% herror = mean(hmax~=hlab)
% 
% hconferr = 0;
% for k = 1:numel(hlab)
%     hconferr = hconferr + 1 - hconf(k, hlab(k));
% end
% hconferr = hconferr / numel(hlab)
% 
% % get same_label error
% %[serror, sconferr] = sameLabelError(spfeatures, spvlabels, sphlabels, ...
% %    {imsegs(:).adjmat}, vclassifier, hclassifier)
% 
% errors.verr = verror;
% errors.vcerr = vconferr;
% errors.herr = herror;
% errors.hcerr = hconferr;
% %errors.serr = serror;
% %errors.scerr = sconferr;
% 
% [hdata2, hlab2] = formatData(spfeatures, spvlabels);
% hconf2 = test_boosted_dt_mc(hclassifier, hdata2);
% hconf2 = 1 ./ (1+exp(-hconf2));
% hconf2 = hconf2 ./ repmat(sum(hconf2, 2), 1, size(hconf2, 2));
% 
% lab = (vlab==1)*1 + (vlab==2)*(hlab2+1) + (vlab==3)*7;
% pEst = [vconf(:, 1) repmat(vconf(:, 2), 1, 5).*hconf2 vconf(:, 3)];
% pTrue = repmat((1:7), numel(lab), 1).*repmat(lab, 1, 7);
% 
% size(pEst)
% size(pTrue)
% 
% [pMean, pVal, estBias] = evaluateProbabilityEstimate(pEst(:), pTrue(:), 0.05)
% pErrors.pMean = pMean;  
% pErrors.pVal = pVal;
% pErrors.estBias = estBias;
% 
% 
% 


% 
% 
% 
% pVal = [step/2:step:(1-step/2)];
% 
% nsteps = numel(pVal);
% 
% estBias = 0;
% 
% for v = 1:numel(pVal)
%     v1 = pVal(v)-step/2;
%     v2 = min(v1+step, 1);
%     
%     ind = find((pEst >= v1) & (pEst <= v2));
%     
%     pMean(v) = mean(pTrue(ind));
%     estBias = estBias + numel(ind)/numel(pEst)*(pMean(v)-pVal(v));
% end