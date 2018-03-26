function [edata, adjlist, perim, boundmap] = mcmcGetEdgeData(imsegs, spdata)
% edata(nadj, nfeatures)
%   edata feature descriptions:
%      01 - 03: abs diff mean rgb
%      04 - 06: abs diff hsv conversion
%      07 - 07: chi-squared hue histogram
%      08 - 08: chi-squared sat histogram
%      09 - 23: abs diff mean texture response
%      24 - 24: chi-squared texture response histogram
%      25 - 26: abs diff mean x-y
%      27 - 27: smaller area / larger area  
%      28 - 28: percent of boundary wrt smaller superpixel
%      29 - 29: straightness of boundary


% spdata(nseg, nf)
%   spdata feature descriptions:
%      01 - 03: mean rgb
%      04 - 06: hsv conversion
%      07 - 11: hue histogram
%      12 - 14: sat histogram
%      15 - 29: mean texture response
%      30 - 44: texture response histogram
%      45 - 46: mean x-y
%      47 - 48: 10th, 90th perc. x
%      49 - 50: 10th, 90th perc. y
%      51 - 51: h / w 
%      52 - 52: area


nfeatures = 29;

nseg = imsegs.nseg;

[boundmap, perim] = mcmcGetSuperpixelBoundaries_fast(imsegs);

nadj = 0;
for s1 = 1:nseg
    nadj = nadj + numel(find(perim(s1, s1+1:end)>0));
end

adjlist = zeros(nadj, 2);
edata = zeros(nadj, nfeatures);
c = 0;
for s1 = 1:nseg
    ns1 = numel(find(perim(s1, s1+1:end)>0));
    adjlist(c+1:c+ns1, 1) = s1;
    adjlist(c+1:c+ns1, 2) = s1 + find(perim(s1, s1+1:end)>0);
    c = c + ns1;
end

[imh, imw] = size(imsegs.segimage);

for k = 1:nadj
    s1 = adjlist(k, 1);
    s2 = adjlist(k, 2);
    
    % differences of spdata
    edata(k, 1:6) = abs(spdata(s1, 1:6) - spdata(s2, 1:6));
    edata(k, 7) = jensenshannon(spdata(s1, 7:11), spdata(s2, 7:11));
    edata(k, 8) = jensenshannon(spdata(s1, 12:14), spdata(s2, 12:14));
    edata(k, 9:23) = abs(spdata(s1, 15:29) - spdata(s2, 15:29));
    edata(k, 24) = jensenshannon(spdata(s1, 30:44), spdata(s2, 30:44));
    edata(k, 25:26) = abs(spdata(s1, 45:46) - spdata(s2, 45:46));
    edata(k, 27) = min(spdata([s1 s2], 52)) / max(spdata([s1 s2], 52));
    
    % boundary information
    edata(k, 28) = perim(s1, s2) / ...
        min(sum(perim(s1, :))+sum(perim(:, s1)), sum(perim(s2, :))+sum(perim(:, s2)));    
    bpix = boundmap{s1, s2};
    bx = double(floor((bpix-1)/imh)+1);
    by = double(mod((bpix-1), imh)+1);    

    edata(k, 29) = sqrt((max(bx)-min(bx)).^2 + (max(by)-min(by)).^2) / perim(s1, s2);
end



function x = chisquare(hist1, hist2)
x = sum((hist1-hist2).^2 ./ hist2);


function d = jensenshannon(hist1, hist2)
% Jensen-Shannon divergence (see wikipedia)
avehist = (hist1 + hist2)/2;
d = 0.5*(sum(hist1.*log(hist1./hist2)) + sum(hist2.*log(hist2./hist1)));