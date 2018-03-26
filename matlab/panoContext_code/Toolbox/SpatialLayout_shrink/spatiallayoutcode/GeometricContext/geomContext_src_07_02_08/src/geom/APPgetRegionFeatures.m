function regionFeatures = APPgetRegionFeatures(image, imsegs, currMap, ...
                                regionNums, spdata, vpdata)
% regionFeatures = APPgetRegionFeatures(image, imsegs, currMap, ...
%                                regionNums, spdata, vpdata)                            
% Computes features for segments given by currMap
%
% Input:
%   image: image to be analyzed
%   imsegs: sp structure
%   currMap: maps sp to regions
%   regionNums: indices of regions to be analyzed
%   spdata: sp data structure from APPgetSpData
%   vpdata: vanishing point data structure
% Output:
%   regionFeatures(nRegions, nFeatures): region features
% 
% Assuming there are 12 oriented filters and 12 texton filters the
% features are as follows:
% 1-12: mean abs oriented filter responses (ofr)
% 13: edginess - mean of all ofrs
% 14: index of largest ofr
% 15: dominance of largest response - max - median of ofr
% textons (would be 16-29) not often used 
%       16-27: mean abs texton filter response (tfr)
%       28: index of largest tfr
%       29: dominance of largest tfr (max - med)
% 16-18: red, green, blue means
% 19-21: hue, saturation, value means (rgb2hsv)
% 22-23: y location, x location means
% 24-29: hue - histogram (5 bins), entropy
% 30-33: sat - histogram (3 bins), entropy
% 34-37: location - 10% and 90% (percentiles) y and x 
% 38-39: edginess center y and x
% 40: number of super-pixels in segmentation
% 41: % of image area region takes up
% 42: number of sides in convex hull polygon
% 43: (num pixels) / (area of convex hull polygon)
% 44: whether the region is contiguous
% 45-60: vanishing point features (see vp2regionfeatures)
% 61-62: 10% and 90% y wrt horizon
% 63   : 0-> below horizon, 1-> straddles horizon, 2-> above horizon
% 64-82: older vanishing point features (see lines_to_vp_info3)
% discrete: 14 44 60 63
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

% WARNING: not using texton features

num_clusters = length(regionNums);

regionFeatures = zeros(num_clusters, 82);

num_angles = size(spdata.orientation,2);
num_textons = 0;

spfeatures = spdata.features; % superpixel features
pixelp = imsegs.npixels ./ sum(imsegs.npixels);
sporient = spdata.orientation;

height = size(image, 1);
width = size(image, 2);

hue = spdata.hue;
sat=  spdata.sat;
yim = 1-repmat([(0:height-1)/(height-1)]', 1, width);
xim = repmat([(0:width-1)/(width-1)], height, 1);

cinds = APPspInds2regionInds(currMap, imsegs.seginds);
edgeim = spdata.edginess;

for c = 1:num_clusters
    
    spind = find(currMap==regionNums(c));
    cpixelp = pixelp(spind);
    %cpixelp = cpixelp / sum(cpixelp); % line added 02/20/06
    
    % features consisting of mean block responses    
    
    % mean orientation filter responses
    rorient = zeros(num_angles, 1);
    for a = 1:num_angles
        regionFeatures(c, a) = sum(spfeatures(spind, a).*cpixelp);
        rorient(a) = sum(sporient(spind, a).*cpixelp);        
    end
    % mean edginess (overall response)
    regionFeatures(c, num_angles+1) = sum(spfeatures(spind, num_angles+1).*cpixelp);
    % most dominant filter response
    [tmp, regionFeatures(c, num_angles+2)] = max(rorient); 
    % dominance of largest filter (max - median)    
    regionFeatures(c, num_angles+3) = max(rorient) - median(rorient);   
	% mean texton response
    for t = 1:num_textons
        regionFeatures(c, num_angles+3+t) = sum(spfeatures(spind, num_angles+t).*cpixelp);
    end
    nf = num_angles+3+num_textons;
    if num_textons>0
        % most dominant texton response
        [tmp, regionFeatures(c, nf+1)] = max(regionFeatures(c, num_angles+3+(1:num_textons))); 
        % dominance of largest texton reponse
        regionFeatures(c, nf+2) = regionFeatures(c, nf+1) - median(regionFeatures(c, num_angles+3+(1:num_textons)));
        nf = nf +2;
    end
        
    % rgb means
    for b = 1:3
        regionFeatures(c, nf+b) = sum(spfeatures(spind, nf+2+b).*cpixelp);
    end
    % hsv
    regionFeatures(c, nf+(4:6)) = rgb2hsv(regionFeatures(c, nf+(1:3)));
    
    % y and x locs
    regionFeatures(c, nf+7) = sum(spfeatures(spind, nf+7).*cpixelp);
    regionFeatures(c, nf+8) = sum(spfeatures(spind, nf+8).*cpixelp);
    
    % get the values that pertain to this region
    rhue = hue(cinds{c});
    rsat = sat(cinds{c});
    rx = xim(cinds{c});
    ry = yim(cinds{c});
    redge = edgeim(cinds{c});
    npix = length(cinds{c});
    
    % hue and sat histograms
    nf = nf+8;
    hue_histogram = (hist(rhue, [0.1:0.2:0.9]) + 0.01)/(length(rhue)+0.05);    
    regionFeatures(c, nf+(1:5)) = hue_histogram / sum(hue_histogram);
    regionFeatures(c, nf+6) = -1*sum(hue_histogram.*log(hue_histogram)) / log(length(hue_histogram));
    sat_histogram = (hist(rsat, [0.167:0.333:0.833]) + 0.01)/(length(rsat)+0.03); 
    regionFeatures(c, nf+(7:9)) = sat_histogram / sum(sat_histogram);
    regionFeatures(c, nf+10) = -1*sum(sat_histogram.*log(sat_histogram)) / log(length(sat_histogram));
    
    % location - 10% and 90% percentiles of y and x
    sorty = sort(ry);
    regionFeatures(c, nf+11) = sorty(ceil(npix/10));
    regionFeatures(c, nf+12) = sorty(ceil(9*npix/10));
    sortx = sort(rx);
    regionFeatures(c, nf+13) = sortx(ceil(npix/10));
    regionFeatures(c, nf+14) = sortx(ceil(9*npix/10));    
    
    % center of edginess y and x
    if regionFeatures(c, nf+11)==regionFeatures(c, nf+12)
        regionFeatures(c, nf+15) = 0.5;
    else
        center_y = sum(redge.*ry) / sum(redge);
        regionFeatures(c, nf+15) = sum((ry < center_y).*ry)/sum(ry);
    end
    if regionFeatures(c, nf+13)==regionFeatures(c, nf+14)
        regionFeatures(c, nf+16) = 0.5;
    else
        center_x = sum(redge.*rx) / sum(redge);
        regionFeatures(c, nf+16) = sum((rx < center_x).*rx)/sum(rx);
    end    

    % num superpixels, % of image area
    regionFeatures(c, nf+17) = length(spind);
    regionFeatures(c, nf+18) = npix / width / height;
    
    % polygon: num sides, area ratio
    try
        [polyi, polya] = convhull(rx, ry);
        regionFeatures(c, nf+19) = length(polyi)-1;
        regionFeatures(c, nf+20) = npix / (polya*width*height);
    catch
        regionFeatures(c, nf+(19:20)) = [4 0.75];
    end

    % whether contiguous
    regionFeatures(c, nf+21) = 1;
    for s = 1:length(spind)
        isadj = imsegs.adjmat(setdiff(spind, spind(s)), spind(s));
        if sum(isadj) == 0
            regionFeatures(c, nf+21) = 0;         
            break;
        end
    end

    % vanishing point features
    nf = nf + 21;
    region_center = [sorty(ceil(npix/2)) sortx(ceil(npix/2))];  

    rbounds = [sortx(1) sortx(end) 1-sorty(end) 1-sorty(1)];
    regionFeatures(c, nf+(1:16)) = ...
        APPvp2regionFeatures(spind, vpdata, region_center, rbounds, imsegs);

    % y-location with respect to estimated horizon
    if ~isnan(vpdata.hpos)
        % location - 10% and 90% percentiles of y and x
        regionFeatures(c, nf+17) = regionFeatures(c, 48) - (1-vpdata.hpos); % bottom 10 pct wrt horizon
        regionFeatures(c, nf+18) = regionFeatures(c, 49) - (1-vpdata.hpos); % top 10 pct wrt horizon
        % 1 -> completely under horizon, 2-> straddles horizon, 3-> completely above horizon
        regionFeatures(c, nf+19) = (regionFeatures(c, nf+20)>0) + (regionFeatures(c, nf+21)>0) + 1;
    else % horizon was not estimated with high confidence
        regionFeatures(c, nf+(17:18)) = regionFeatures(c, [48:49])-0.5;
        regionFeatures(c, nf+19) = 4;  % signifies no data-estimated horizon
    end
    
    [h w] = size(image);
    region_center(1) = region_center(1) - 0.5;
    region_center(2) = region_center(2) - 0.5;    
    regionFeatures(c, nf+(20:38)) = ...
        APPgetVpFeatures(vpdata.spinfo(spind), vpdata.lines, region_center, [h w]);    
      
end
