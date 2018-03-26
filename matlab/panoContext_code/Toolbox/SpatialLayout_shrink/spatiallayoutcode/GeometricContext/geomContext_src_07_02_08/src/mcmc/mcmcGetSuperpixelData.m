function features = mcmcGetSuperpixelData(im, imsegs)
% features(nseg, nf)
%   Feature Descriptions:
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

ncolor = 6; % mean rgb, hsv
ncolhist = 8;
ntexthist = 15;
nloc = 8; % mean x-y, 10th, 90th percentile x-y, w/h, area
filtext = makeLMfilters;
ntext = size(filtext, 3);
nseg = imsegs.nseg;
npix = imsegs.imsize(1)*imsegs.imsize(2);
if ntext~=15
    disp('warning: expecting 15 texture filters')
end

features = zeros(nseg, ncolor+ntext+nloc+ncolhist+ntexthist);

% color
imhsv = rgb2hsv(im);
grayim = imhsv(:, :, 3);
[imh imw] = size(grayim);
hhist = 1 + (imhsv(:, :, 1) > 0.2) + (imhsv(:, :, 1) > 0.4) + ...
    (imhsv(:, :, 1) > 0.6) + (imhsv(:, :, 1) > 0.8);
shist = 1 + (imhsv(:, :, 2) > 0.33) + (imhsv(:, :, 2) > 0.67);

% texture
imtext = zeros([imh imw ntext]);
for f = 1:ntext
    imtext(:, :, f) = abs(imfilter(im2single(grayim), filtext(:, :, f), 'same'));    
end
[tmp, texthist] = max(imtext, [], 3);

% location
yim = 1-repmat(((0:imh-1)/(imh-1))', 1, imw);
xim = repmat(((0:imw-1)/(imw-1)), imh, 1); 

spstats = regionprops(imsegs.segimage, 'PixelIdxList');

for s = 1:nseg    
    % rgb means
    f = 0;
    for k = 1:3
        features(s, f+k) = mean(im(spstats(s).PixelIdxList+(k-1)*imw*imh));
    end
    f = f + 3;
    % hsv means
    features(s, f+[1:3]) = rgb2hsv(features(s, [1:3]));
    f = f + 3;
    % hue histogram    
    features(s, f+[1:5]) = hist(hhist(spstats(s).PixelIdxList), [1:5])+1;
    features(s, f+[1:5]) = features(s, f+[1:5]) / sum(features(s, f+[1:5]));
    f = f + 5;
    % sat histogram        
    features(s, f+[1:3]) = hist(shist(spstats(s).PixelIdxList), [1:3])+1;
    features(s, f+[1:3]) = features(s, f+[1:3]) / sum(features(s, f+[1:3]));
    f = f + 3;
    % texture means
    for k = 1:ntext
        features(s, f+k) = mean(imtext(spstats(s).PixelIdxList+(k-1)*npix));
    end
    f = f + ntext;
    % texture histogram
    features(s, f+(1:ntext)) = hist(texthist(spstats(s).PixelIdxList), (1:ntext))+1;
    features(s, f+(1:ntext)) = features(s, f+(1:ntext)) / sum(features(s, f+(1:ntext)));
    f = f + ntext;
    % location
    xvals = xim(spstats(s).PixelIdxList);
    yvals = yim(spstats(s).PixelIdxList);
    features(s, f+1) = mean(xvals);
    features(s, f+2) = mean(yvals);
    sxvals = sort(xvals);
    syvals = sort(yvals);
    features(s, f+3) = sxvals(ceil(numel(sxvals)/10));
    features(s, f+4) = sxvals(ceil(9*numel(sxvals)/10));
    features(s, f+5) = syvals(ceil(numel(syvals)/10));
    features(s, f+6) = syvals(ceil(9*numel(syvals)/10));    
    features(s, f+7) = (features(s, f+4) - features(s, f+3)) / ...
        (features(s, f+6) - features(s, f+5) + 0.0001);
    features(s, f+8) = imsegs.npixels(s)/imh/imw;
end

