function [errors, perrors] = ...
    mcmcTestSuperpixelClassifier(spfeatures, imsegs, vclassifier, hclassifier)

spvlabels = {imsegs(:).vert_labels};
sphlabels = {imsegs(:).horz_labels};

% test vertical
[vdata, vlab, imind] = formatData(spfeatures, spvlabels);
vconf = test_boosted_dt_mc(vclassifier, vdata);
vconf = 1 ./ (1+exp(-vconf));
vconf = vconf ./ repmat(sum(vconf, 2), 1, size(vconf, 2));
[tmp, vmax] = max(vconf, [], 2);
verror = mean(vmax~=vlab)

vconferr = 0;
for k = 1:numel(vlab)
    vconferr = vconferr + 1 - vconf(k, vlab(k));
end
vconferr = vconferr / numel(vlab)

% test subclass
[hdata, hlab] = formatData(spfeatures, sphlabels);
hconf = test_boosted_dt_mc(hclassifier, hdata);
hconf = 1 ./ (1+exp(-hconf));
hconf = hconf ./ repmat(sum(hconf, 2), 1, size(hconf, 2));
[tmp, hmax] = max(hconf, [], 2);
herror = mean(hmax~=hlab)

hconferr = 0;
for k = 1:numel(hlab)
    hconferr = hconferr + 1 - hconf(k, hlab(k));
end
hconferr = hconferr / numel(hlab)

% get same_label error
%[serror, sconferr] = sameLabelError(spfeatures, spvlabels, sphlabels, ...
%    {imsegs(:).adjmat}, vclassifier, hclassifier)

errors.verr = verror;
errors.vcerr = vconferr;
errors.herr = herror;
errors.hcerr = hconferr;
%errors.serr = serror;
%errors.scerr = sconferr;

[hdata2, hlab2, imind] = formatData(spfeatures, spvlabels);
hconf2 = test_boosted_dt_mc(hclassifier, hdata2);
hconf2 = 1 ./ (1+exp(-hconf2));
hconf2 = hconf2 ./ repmat(sum(hconf2, 2), 1, size(hconf2, 2));

lab = (vlab==1)*1 + (vlab==2).*(hlab2+1) + (vlab==3)*7;

pEst = [vconf(:, 1) repmat(vconf(:, 2), 1, 5).*hconf2 vconf(:, 3)];
pTrue = (repmat((1:7), numel(lab), 1)==repmat(lab, 1, 7));

[pMean, pVal, estBias] = evaluateProbabilityEstimate(pEst(:), pTrue(:), 0.05);
perrors.pMean = pMean;  
perrors.pVal = pVal;
perrors.estBias = estBias;

perrors.pMean2 = 0;
perrors.estBias2 = 0;
snv = 0;
for k = 1:5
    testind = find((imind > (1-k)*50) & (imind < k*50));
    trainind = setdiff((1:numel(imind)),testind)';        
    
    [pMean, pVal, estBias] = evaluateProbabilityEstimate(pEst(trainind), pTrue(trainind), 0.05);
    perrors.pMeanTr{k} = pMean;
    
    pEstTest = pEst(testind);
    for i = 1:numel(pEstTest)
        [val, maxind] = min(abs(pVal-pEstTest(i)));
        pEstTest(i) = pMean(maxind);
    end
%    disp(num2str(numel(pEstTest)))
    [pMean, pVal, estBias, nv] = evaluateProbabilityEstimate(pEstTest, pTrue(testind), 0.05);
    snv = snv + nv;    
    perrors.pMean2 = perrors.pMean2 + pMean.*nv;
    perrors.estBias2 = perrors.estBias2 + estBias/5;
end
perrors.pMean2 = perrors.pMean2 ./ max(snv, 1E-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data, lab, imind] = formatData(features, labels, labels2)
% concatenate data and labels, and give adjacency information

nimages = numel(features);
[tmp, nvars] = size(features{1});
% count superpixels
nsp = 0;
for f = 1:nimages
    nsp = nsp + sum(labels{f}~=0);
end
data = zeros(nsp, nvars);
lab = zeros(nsp, 1);
imind = zeros(nsp, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = find(labels{f}~=0);
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    lab(vc+1:vc+numel(ind)) = labels{f}(ind);
    imind(vc+1:vc+numel(ind)) = f;
    vc = vc + numel(ind);    
end



function [serror, sconferr] = sameLabelError(features, vlabels, hlabels, adjmat, vclassifier, hclassifier)

nimages = numel(features);
[tmp, nvars] = size(features{1});
% count superpixels
nsp = 0;
for f = 1:nimages
    nsp = nsp + sum((vlabels{f}~=0) & ((hlabels{f}~=0) | (vlabels{f}~=2)));
end
data = zeros(nsp, nvars);
lab = zeros(nsp, 1);
adj = cell(nsp, 1);
issame = cell(nsp, 1);

% concatenate data
vc = 0;
nadj = 0;
for f = 1:nimages
    ind = find((vlabels{f}~=0) & ((hlabels{f}~=0) | (vlabels{f}~=2)));
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    for k = 1:numel(ind)
        restind = ind(k+1:end);
        adj{vc+k} = restind(find(adjmat{f}(ind(k), restind)));
        issame{vc+k} = zeros(numel(adj{vc+k}), 1);
        for k2 = 1:numel(adj{vc+k})
            issame{vc+k}(k2) = (vlabels{f}(ind(k))==vlabels{f}(adj{vc+k}(k2))) && ...
                (hlabels{f}(ind(k))==hlabels{f}(adj{vc+k}(k2)));
        end
        nadj = nadj  +numel(adj{vc+k});
    end
    vc = vc + numel(ind);
end
            
vconf = test_boosted_dt_mc(vclassifier, data);    
vconf = 1 ./ (1+exp(-vconf));
vconf = vconf ./ repmat(sum(vconf, 2), 1, size(vconf, 2));
hconf = test_boosted_dt_mc(hclassifier, data);
hconf = 1 ./ (1+exp(-hconf));
hconf = hconf ./ repmat(sum(hconf, 2), 1, size(hconf, 2));

sconf = zeros(nadj, 1);
slab = zeros(nadj, 1);
sc = 0;
for k = 1:numel(adj)
    for k2 = 1:numel(adj{k})
        sc = sc + 1;
        slab(sc) = issame{k}(k2);
        sconf(sc) = sum(vconf(k, [1 3]).*vconf(k2, [1 3])) + ...
            sum((vconf(k, 2)*hconf(k, :)) .* (vconf(k2, 2)*hconf(k2, :)));
    end
end

ind1 = find(slab==1);
ind2 = find(slab==0);
disp(num2str([mean(sconf(ind1)) mean(sconf(ind2))]))
[sorted] = sort(sconf);
tval = sorted(floor(numel(sconf)*(1-mean(slab))));
serror = mean(slab~=(sconf>tval));
sconferr = mean(abs(slab-sconf));

