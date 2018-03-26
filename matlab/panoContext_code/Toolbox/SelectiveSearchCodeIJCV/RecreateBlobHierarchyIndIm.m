function [hBlobs blobsInit blobsRest] = RecreateBlobHierarchyIndIm(blobIndIm, blobBoxes, hierarchy)
% function hBlobs = RecreateBlobHierarchyIndIm(blobIndIm, boxes, hierarchy)
%
% Recreate hierarchy from the initial segmentation image
%
% blobIndIm:            Image with indices denoting segments
% blobBoxes:            Boxes belonging to blobs in blobIndIm
% hierarchy:            Hierarchy denoting hierarchical merging
% 
% hBlobs:               All blobs in the hierarchy
% blobsInit:            The initial blobs
% blobsRest:            All blobs but the initial blobs
%
%     Jasper Uijlings - 2013

% Get blobs of initial segmentation
blobsInit = SegmentIndices2Blobs(blobIndIm, blobBoxes);

% Add sizes
blobsInit = BlobAddSizes(blobsInit);

% Reconstruct hierarchy
hBlobs = RecreateBlobHierarchy(blobsInit, hierarchy);

if nargout == 3
    blobsRest = hBlobs(length(blobsInit)+1:end);
end