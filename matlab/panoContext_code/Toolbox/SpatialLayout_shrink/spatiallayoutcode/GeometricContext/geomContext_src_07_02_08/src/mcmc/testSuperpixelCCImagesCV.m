function [vacc, hacc, vcm, hcm] = testSuperpixelCCImagesCV(imsegs, vclassifier, hclassifier, segdata, smap, ncv)

vacc = 0;
hacc = 0;
vtotal = 0;
htotal = 0;

vcm = zeros(3);
hcm = zeros(5);

nimages = numel(imsegs);

for f = 1:nimages  
    
    k = ceil(f / (nimages/ncv));
    
    [pv, ph] = segClassify(segdata{f}, vclassifier(k), hclassifier(k), smap{f});

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

vcm = vcm ./ repmat(sum(vcm, 2), 1, size(vcm, 2));
hcm = hcm ./ repmat(sum(hcm, 2), 1, size(hcm, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [pv, ph] = segClassify(segdata, vclassifier, hclassifier, smap)

nsp = numel(smap);
nseg = max(smap);

pv = zeros(nsp, 3);
ph = zeros(nsp, 5);

for s = 1:nseg

    ind = find(smap==s);

    % probability of superpixel main labels
    tpv = test_boosted_dt_mc(vclassifier, segdata(s, :));
    tpv = 1 ./ (1+exp(-tpv));
    tpv = tpv ./ repmat(sum(tpv, 2), 1, size(tpv, 2));
    pv(ind, :) = repmat(tpv, numel(ind), 1);    
    
    % probability of superpixel sub labels
    tph = test_boosted_dt_mc(hclassifier, segdata(s, :));
    tph = 1 ./ (1+exp(-tph));
    tph = tph ./ repmat(sum(tph, 2), 1, size(tph, 2));
    ph(ind, :) = repmat(tph, numel(ind), 1);
end