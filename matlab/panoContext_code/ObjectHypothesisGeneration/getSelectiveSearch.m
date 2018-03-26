function [ SS1 ] = getSelectiveSearch( rotImg )
%GETSELECTIVESEARCH Get hierarchy grouping of segmentation under different
%color space.
%   Detailed explanation goes here

img1 = rotImg;

%% selective search, object proposal generation
% color space
colorTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};
% Here you specify which similarity functions to use in merging
simFunctionHandles = {@SSSimColourTextureSizeFillOrig, @SSSimTextureSizeFill, @SSSimBoxFillOrig, @SSSimSize};
simFunctionHandles = simFunctionHandles(1:2); % Two different merging strategies

% Thresholds for the Felzenszwalb and Huttenlocher segmentation algorithm.
% Note that by default, we set minSize = k, and sigma = 0.8.
k = 200; % controls size of segments of initial segmentation. 
minSize = k;
sigma = 0.8;

SS1 = repmat(struct('boxes',[],'blobIndIm',[],'blobBoxes',[],'hierarchy',[]), length(colorTypes), 1);
for i = 1:length(colorTypes)
    colorType = colorTypes{i};
    [SS1(i).boxes, SS1(i).blobIndIm, SS1(i).blobBoxes, SS1(i).hierarchy] = ...
        Image2HierarchicalGrouping(im2uint8(img1), sigma, k, minSize, colorType, simFunctionHandles);
end


end

