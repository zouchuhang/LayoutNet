function [abo mabo boScores avgNumSegments] = BlobAverageBestOverlap(gtBlobs, gtNrs, blobIndIm, blobBoxes, hierarchy, minWidth)
% [abo mabo boScores avgNumSegments] = BlobAverageBestOverlap(gtBlobs,
%                   gtNrs, blobIndIm, blobBoxes, hierarchy, minWidth)
%
% Calculate Average Best Overlap scores
%
% gtBlobs:      Cell array of ground truth segments per class (see
%               GetAllObjectBoxes)
% gtNrs:        Cell array with image nrs corresponding to ground truth.
% blobIndIm:    Image with indices per blob (mexFelzenSegmentIndex)
% blobBoxes:    Boxes corresponding to blobs in blobIndIm
% hierarchy:    Hierarchy necessary to reconstruct all blobs in grouping
% minWidth:     (optional) Filter out blobs with a width smaller than minWidth.
%
% abo:          Average Best Overlap per class (Pascal Overlap criterion)
% mabo:         Mean Average Best Overlap (mean(abo))
% boScores:     Best Overlap Score per GT segment.
% avgNumBlobs:  Average number of blobs per image
%
%     Jasper Uijlings - 2013

if ~exist('minWidth', 'var')
    minWidth = 0;
end

nClasses = length(gtBlobs);

% Memory initialization
numSegments = zeros(length(blobIndIm), 1);
boScores = cell(1, nClasses);
for cI = 1:nClasses
    boScores{cI} = length(gtBlobs{cI});
end

% indices per class
classIdx = ones(1, nClasses);

for cI=1:length(gtBlobs)
    for i=1:length(gtBlobs{cI})
        testImNr = gtNrs{cI}(i);
        
        % the hierarchy here contains possibly multiple groupings with
        % different initial measures
        testBlobsT = cell(length(hierarchy{testImNr}), 1);
        testBlobsT{1} = RecreateBlobHierarchyIndIm(blobIndIm{testImNr}, blobBoxes{testImNr}, hierarchy{testImNr}{1});
        for j=2:length(hierarchy{testImNr}) % Without initial blobs here
            [aa bb testBlobsT{j}] = RecreateBlobHierarchyIndIm(blobIndIm{testImNr}, blobBoxes{testImNr}, hierarchy{testImNr}{j});
        end
        testBlobs = cat(1, testBlobsT{:});
        
        % Get rid of too small blobs
        testBlobs = FilterBlobsWidth(testBlobs, minWidth);
        numSegments(testImNr) = length(testBlobs);        
      
        % Calculate overlap scores
        boScores{cI}(classIdx(cI)) = BlobBestOverlap(testBlobs, gtBlobs{cI}(i));
        
        classIdx(cI) = classIdx(cI) + 1;
    end
end

abo = zeros(nClasses, 1);

for cI = 1:nClasses
    abo(cI) = mean(boScores{cI});
end

mabo = mean(abo);    

% Average of numSegments. Make sure that only images for which the
% numSegments are actually calculated are taken into account.
avgNumSegments = mean(numSegments(numSegments > 0));
    