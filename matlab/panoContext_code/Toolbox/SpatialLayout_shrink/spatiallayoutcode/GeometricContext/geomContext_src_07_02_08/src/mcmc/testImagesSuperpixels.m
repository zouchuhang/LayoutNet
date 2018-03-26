function [vacc, hacc, vcm, hcm] = testImagesSuperpixels(imsegs, vclassifierSP, hclassifierSP, spdata)

vacc = 0;
hacc = 0;
vtotal = 0;
htotal = 0;

vcm = zeros(3);
hcm = zeros(5);

for f = 1:numel(imsegs)  
    [pv, ph] = spClassify(spdata{f}, vclassifierSP, hclassifierSP);

    [vmaxval, vmax] = max(pv, [], 2);
    [hmaxval, hmax] = max(ph, [], 2);   
    
    vlab = imsegs(f).vert_labels(:);
    hlab = imsegs(f).horz_labels(:);
    npixels = imsegs(f).npixels;
    npixels = npixels/sum(npixels);

    vacc = vacc + sum((vmax==vlab).*npixels);
    vtotal = vtotal + sum((vlab~=0).*npixels);
    hacc = hacc + sum((hmax==hlab).*npixels);
    htotal = htotal + sum((hlab~=0).*npixels);   

    for s = 1:numel(vmax)
        if vlab(s)~=0
            vcm(vlab(s), vmax(s)) = vcm(vlab(s), vmax(s)) + npixels(s)/numel(imsegs);
        end
        if hlab(s)~=0
            hcm(hlab(s), hmax(s)) = hcm(hlab(s), hmax(s)) + npixels(s)/numel(imsegs);
        end        
    end
end

vacc = vacc / vtotal;
hacc = hacc / htotal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [pvSP, phSP] = spClassify(spdata, vclassifierSP, hclassifierSP)

% probability of superpixel main labels
pvSP = test_boosted_dt_mc(vclassifierSP, spdata);
pvSP = 1 ./ (1+exp(-pvSP));
pvSP = pvSP ./ repmat(sum(pvSP, 2), 1, size(pvSP, 2));

% probability of superpixel sub labels
phSP = test_boosted_dt_mc(hclassifierSP, spdata);
phSP = 1 ./ (1+exp(-phSP));
phSP = phSP ./ repmat(sum(phSP, 2), 1, size(phSP, 2));