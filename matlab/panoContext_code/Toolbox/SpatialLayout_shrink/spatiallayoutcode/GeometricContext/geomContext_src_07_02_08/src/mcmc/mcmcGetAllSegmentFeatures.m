function features = mcmcGetAllSegmentFeatures(imsegs, imdir, gtmaps, spfeatures, vclassifierSP, hclassifierSP)

for f = 1:numel(imsegs)       
    disp([num2str(f) ': ' imsegs(f).imname])
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    %[pvSP, phSP] = spClassify(spfeatures{f}, vclassifierSP, hclassifierSP);
    imdata = mcmcComputeImageData(im, imsegs(f));
    %imdata.pvSP = pvSP;
    %imdata.phSP = phSP;        
    features{f} = mcmcGetSegmentFeatures(imsegs(f), spfeatures{f}, imdata, gtmaps{f}, [1:max(gtmaps{f})]);
end



function [pvSP, phSP] = spClassify(spdata, vclassifierSP, hclassifierSP)

% probability of superpixel main labels
pvSP = test_boosted_dt_mc(vclassifierSP, spdata);
pvSP = 1 ./ (1+exp(-pvSP));
pvSP = pvSP ./ repmat(sum(pvSP, 2), 1, size(pvSP, 2));

% probability of superpixel sub labels
phSP = test_boosted_dt_mc(hclassifierSP, spdata);
phSP = 1 ./ (1+exp(-phSP));
phSP = phSP ./ repmat(sum(phSP, 2), 1, size(phSP, 2));
