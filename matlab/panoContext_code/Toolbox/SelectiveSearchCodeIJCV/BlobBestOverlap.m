function [scores index] = BlobBestOverlap(gtBlobs, testBlobs)
% [scores index] = BlobBestOverlap(gtBlobs, testBlobs)
% 
% Get overlap scores (Pascal-wise) for test  blobs
%
% groundTruthBlob:         ground truth blobs
% test:                    Test blobs
%
% scores:                  Highest overlap scores for each test blob.
% index:                   Index for each test blob which ground truth blob
%                          is best
%
%     Jasper Uijlings - 2013

numTarget = length(gtBlobs);
numTest = length(testBlobs);

scoreM = zeros(numTest, numTarget);

for i=1:numTest
    for j=1:numTarget
        scoreM(i,j) = PascalOverlapBlob(gtBlobs{j}, testBlobs{i});
    end
end

[scores index] = max(scoreM, [], 2);

