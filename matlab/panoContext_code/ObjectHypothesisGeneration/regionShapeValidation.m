function [ hBlobs_valid, hBlobs_centroid ] = regionShapeValidation( hBlobs, H, W, min_size, aspectratio, arearatio, iscuboid )
%REGIONSHAPEVALIDATION Summary of this function goes here
%   check whether to find a cuboid or rectangle fitting based on
%   segmentation shape
numSeg = length(hBlobs);
invalidMask = false(H,W);
if iscuboid
    invalidMask(:,1) = true;   invalidMask(:,W) = true;
    invalidMask(1,:) = true;   invalidMask(H,:) = true;
    invalidMask(round(H/2), 1) = true;
    invalidMask(round(H/2), round(W/4)) = true;
    invalidMask(round(H/2), round(W/2)) = true;
    invalidMask(round(H/2), round(3*W/4)) = true;
    invalidMask(round(H/2), round(W)) = true;
    invalidMask = bwdist(invalidMask)<10;
end
SE = strel('disk', 5, 4);
SS = strel('disk', 2, 4);

hBlobs_valid = false(numSeg, 1);
hBlobs_centroid = zeros(numSeg,2);

parfor bid = 1:length(hBlobs)
%     fprintf('%d/%d\n', bid, length(hBlobs));
    
    if hBlobs{bid}.size<min_size
        continue;
    end
    seg = false(H,W);
    rect = hBlobs{bid}.rect;
    seg(rect(1):rect(3),rect(2):rect(4)) = hBlobs{bid}.mask;
    
    binaryMap = imdilate( seg, SE);
    binaryMap = imerode( binaryMap, SE);
    binaryMap = imerode( binaryMap, SS);
    binaryMap = imdilate( binaryMap, SS);
    binaryMap = imfill(binaryMap, 'hole');
    
    CC = bwconncomp(binaryMap);
    valid = true(CC.NumObjects,1);
    for j = 1:CC.NumObjects
        if length(CC.PixelIdxList{j})<min_size
            valid(j) = false;
        end
        if any(invalidMask(CC.PixelIdxList{j}))
            valid(j) = false;
        end
    end
    CC.NumObjects = sum(valid);
    CC.PixelIdxList = CC.PixelIdxList(valid);
    
    if CC.NumObjects ~= 1
        continue;
    end
    
    R = regionprops( CC, 'BoundingBox', 'Area', 'ConvexArea', 'Centroid');
    boundingbox = R(1).BoundingBox;
    asrt = boundingbox(4)/boundingbox(3);
    if asrt>aspectratio || asrt<1/aspectratio
        continue;
    end
    if  R(1).Area/(R(1).ConvexArea+0.0001)<arearatio
        continue;
    end
    hBlobs_valid(bid) = true;
    hBlobs_centroid(bid,:) = R(1).Centroid;
end



end

