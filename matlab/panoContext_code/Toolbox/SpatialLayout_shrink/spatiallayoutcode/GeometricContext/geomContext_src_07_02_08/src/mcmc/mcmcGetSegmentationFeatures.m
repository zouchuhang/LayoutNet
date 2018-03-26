function data = mcmcGetSegmentationFeatures(pvSP, phSP, pE, adjlist, npixels, smap, sind)
% features:
% 01 - 07: sorted mean superpixel label probabilities (main+subclass)
% 08 - 10: sorted mean superpixel main label probabilities 
% 11 - 17: sorted max superpixel label probabilities (main+subclass)
% 18 - 20: sorted max superpixel main label probabilities 
% 21 - 27: sorted histogram of most likely labels (main+subclass)
% 28 - 30: sorted histogram of most likely main labels 
% 31 - 32: entropy of labels / main labels (according to initial est)
% 33 - 34: mean and min of pairwise same-label likelihoods
% 35: number of superpixels
% 36: area

nseg = numel(sind);

data = zeros(nseg, 36);

totalpixels = sum(npixels);

c = 0;

for s = sind(:)'
    
    c = c + 1;
    
    spind = find(smap==s);
    
    if numel(spind)==0
        disp([s max(smap)])
        disp('error: empty segment')
    end
    
    pv = pvSP(spind, :);
    ph = phSP(spind, :);
    
    plab = [pv(:, 1) repmat(pv(:, 2), [1 size(ph, 2)]).*ph pv(:, 3)];
 
    % weighted (by npixels) mean of confidences
    npix = npixels(spind)/sum(npixels(spind));
    for k = 1:3
        meanpv(k) = sum(pv(:, k).*npix);
    end
    for k = 1:7
        meanplab(k) = sum(plab(:, k).*npix);
    end       
    
    %meanpv = mean(pv, 1);
    %meanplab = mean(plab, 1);
    
    % highest confidence for each label
    [mlvVal, tmp] = max(pv, [], 1);
    [mllabVal, tmp] = max(plab, [], 1);
    
    data(c, 1:7) = sort(meanplab, 'descend');
    data(c, 8:10) = sort(meanpv, 'descend');
    
    data(c, 11:17) = sort(mllabVal, 'descend');
    data(c, 18:20) = sort(mlvVal, 'descend');
 
    % most likely label for each superpixel
    [tmp, mlv] = max(pv, [], 2);
    [tmp, mllab] = max(plab, [], 2);
    histv = hist(mlv, [1:3]);
    histlab = hist(mllab, [1:7]);
    data(c, 21:27) = sort(histlab, 'descend');
    data(c, 28:30) = sort(histv, 'descend');    
        
    histv = histv / sum(histv);        
    histlab = histlab / sum(histlab);
    
    % entropy avoiding divide by zero
    ind = find(histlab>0);
    data(c, 31) = -sum(log(histlab(ind)).*histlab(ind))/log(2);
    ind = find(histv>0);
    data(c, 32) = -sum(log(histv(ind)).*histv(ind))/log(2);

    
    ind = find((smap(adjlist(:, 1))==s) & (smap(adjlist(:, 2))==s));    
    pEs = pE(ind); 
    if ~isempty(pEs)
        data(c, 33) = mean(pEs);
        data(c, 34) = min(pEs);
    else
        data(c, 33:34) = 1;
    end
    
    data(c, 35) = numel(spind);
    
    data(c, 36) = sum(npixels(spind))/totalpixels;
end
    
            
    

