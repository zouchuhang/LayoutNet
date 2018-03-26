function spdata = APPgetSpData(image, doogFilters, textonFilters, imsegs)
% spdata = APPgetSpData(image, doogFilters, textonFilters, imsegs)
%
% Gets the features and data corresponding to the superpixels given by
% imsegs.
% features:
%   1-12: mean absolute filter response of diff of oriented Gaussian filters
%   13: mean of 1-12
%   14: argmax of 1-12
%   15: max-median of 1-12
%   possibly texton features (15+1:15+nTextons+2)
%   16-18: mean rgb values
%   19-21: hsv conversion from mean rgb values
%   22-23: mean x and y locations
%
% Input:
%   image: rgb image to be analyzed
%   doogFilters: filters for texture (could be empty)
%   textonFilters: filters for texture (could be empty)
%   imsegs: superpixel structure
% Output:
%   spdata: structure for data corresponding to each superpixel
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005
      
nAngles = size(doogFilters, 3);
        
% compute texture responses
grayim = rgb2gray(image);
    
%disp(num2str(size(image, 1)*size(image, 2)*nAngles /1024/1024*8))
orientationImages = zeros(size(image, 1), size(image, 2), nAngles);
% get and record the filter responses for each angle
for (i=1:nAngles)        
    orientationImages(:, :, i) = abs(conv2(grayim,doogFilters(:, :, i),'same')-grayim);        
end
%spdata.orientationImages = orientationImages;

textonImages = compute_texton_response(grayim, textonFilters);
nTextons = numel(textonFilters);

[height, width, nb] = size(image);

% features: angle response peaks, mean angle response, angle with 
% largest response, dominance of largest filter, mean abs texton
% repsonses, largest texton response, dominance of largest response, 
% mean of rgb hsv x y
nfeatures = nAngles + 3 + nTextons + 6 + 2;
if nTextons > 0
    nfeatures = nfeatures + 2;
end

% for each segmentation
for i = 1:length(imsegs)
    
    nseg = imsegs(i).nseg;
    spdata(i).npixels = imsegs(i).npixels;
    spdata(i).adjmat = imsegs(i).adjmat;
    
    features = zeros(nseg, nfeatures);
               
    % oriented filters
    maar = zeros(nseg, nAngles);
    for a = 1:nAngles
        maar(:, a) = APPgetSpMeans(imsegs(i), orientationImages(:, :, a));
    end
    features(:, 1:nAngles) = maar;
    features(:, nAngles+1) = mean(maar, 2);
    [maxval, features(:, nAngles+2)] = max(maar, [], 2);
    features(:, nAngles+3) = maxval - median(maar, 2);    
    spdata(i).orientation = maar;
    spdata(i).edginess = mean(orientationImages, 3);
    clear orientationImages;
    
    % textons    
    cf = nAngles+3;
    for t = 1:nTextons
        features(:, cf+t) = APPgetSpMeans(imsegs(i), textonImages(:, :, t));
    end
    cf = cf + nTextons;
    if nTextons>0
        [maxval, features(:, cf+1)] = max(features(:, cf+1:cf+nTextons), [], 2);
        features(:, cf+2) = maxval - median(features(:, cf+1:cf+nTextons), 2); 
        clear textonImages;
        cf = cf + 2;
    end            
    
    % color
    rgb = zeros(nseg, 3);
    for b = 1:3
        rgb(:, b) = APPgetSpMeans(imsegs(i), image(:, :, b));
    end
    features(:, cf + (1:3)) = rgb;
    features(:, cf + (4:6)) = rgb2hsv(rgb);
    
    % location
    cf = cf + 6;
    yim = 1-repmat(((0:height-1)/(height-1))', 1, width);
    xim = repmat(((0:width-1)/(width-1)), height, 1);
    features(:, cf+1) = APPgetSpMeans(imsegs(i), yim);
    features(:, cf+2) = APPgetSpMeans(imsegs(i), xim);
    spdata(i).features = features;    
end
    
    
