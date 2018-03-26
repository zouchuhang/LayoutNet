function [ viewPoint, labelID, segCentroid] = getValidViewPoint( segment, minSize )
%GETVALIDVIEWPOINT Summary of this function goes here
%   Detailed explanation goes here
% segment = segment + 1;
numSeg = max(segment(:));
n = histc(segment(:),1:numSeg);
[H,W] = size(segment);

invalidMask = false(H,W);
invalidMask(:,1) = true;   invalidMask(:,W) = true;
invalidMask(1,:) = true;   invalidMask(H,:) = true;
invalidMask(round(H/2), 1) = true;
invalidMask(round(H/2), round(W/4)) = true;
invalidMask(round(H/2), round(W/2)) = true;
invalidMask(round(H/2), round(3*W/4)) = true;
invalidMask(round(H/2), round(W)) = true;
invalidMask = bwdist(invalidMask)<10;


% boundSeg = segment;
% boundSeg(2:H-1, 2:W-1) = 0;
% m = histc(segment(invalidMask(:)), 1:numSeg);
% boundID = find(m);

segValid = n>minSize;
% segValid(m>0) = false;
counter = 0;
segCentroid = zeros(1000,2);
labelID = zeros(1000,1);
% segment = imresize(segment, 0.5);
% segment = segment>0.5;
SE = strel('disk', 5, 4);
SS = strel('disk', 2, 4);

% conditionMask = zeros(H,W);
% I = ismember(segment, find(~segValid));
% conditionMask(I) = 1;

for i = find(segValid)'
    seg = segment == i;
    
    binaryMap = imdilate( seg, SE);
    binaryMap = imerode( binaryMap, SE);
    binaryMap = imerode( binaryMap, SS);
    binaryMap = imdilate( binaryMap, SS);
    binaryMap = imfill(binaryMap, 'hole');
    
    CC = bwconncomp(binaryMap);
    valid = true(CC.NumObjects,1);
    for j = 1:CC.NumObjects
        if length(CC.PixelIdxList{j})<minSize
            valid(j) = false;
%             conditionMask(CC.PixelIdxList{j}) = 1; % too small
        end
        if any(invalidMask(CC.PixelIdxList{j}))
            valid(j) = false;
%             conditionMask(CC.PixelIdxList{j}) = 2; % hit invalid region
        end
    end
    CC.NumObjects = sum(valid);
    CC.PixelIdxList = CC.PixelIdxList(valid);
    
    R = regionprops( CC, 'BoundingBox', 'Area', 'ConvexArea', 'Centroid');
    for j = 1:length(R)
        boundingbox = R(j).BoundingBox;
        asrt = boundingbox(4)/boundingbox(3);
        if asrt>3 || asrt<0.333
%             conditionMask(CC.PixelIdxList{j}) = 3; % aspect ratio
            continue;
        end
        if  R(j).Area/(R(j).ConvexArea+0.0001)>0.65
            counter = counter + 1;
            segCentroid(counter,:) = R(j).Centroid;
            labelID(counter) = i;
        else
%             conditionMask(CC.PixelIdxList{j}) = 4; % convex
        end
    end
    
end
segCentroid = segCentroid(1:counter,:);
labelID = labelID(1:counter);
viewPoint = coords2uv( segCentroid, W, H );

% figure; imshow(segment,[]); hold on
% for i = 1:size(segCentroid,1)
%          plot(segCentroid(i,1), segCentroid(i,2), '--rs','LineWidth',2,...
%                  'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','b',...
%                  'MarkerSize',10);
% end
end

