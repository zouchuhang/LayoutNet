function [abo mabo boScores avgNumBoxes] = BoxAverageBestOverlap(gtBoxes, gtNrs, testBoxes)
% [abo mabo boScores avgNumBoxes] = BoxAverageBestOverlap(gtBoxes, gtNrs, testBoxes)
%
% Calculate Average Best Overlap scores
%
% gtBoxes:      Cell array of ground truth boxes per class (see
%               GetAllObjectBoxes)
% gtNrs:        Cell array with image nrs corresponding to ground truth.
% testBoxes:    Cell array of testboxes per image.
%
% abo:          Average Best Overlap per class (Pascal Overlap criterion)
% mabo:         Mean Average Best Overlap (mean(abo))
% boScores:     Best Overlap Score per GT box.
% avgNumBoxes:  Average number of boxes per image
%
%     Jasper Uijlings - 2013

% Check nr of gt elements
nClasses = length(gtBoxes);

boScores = cell(1, nClasses);
for cI = 1:nClasses
    boScores{cI} = zeros(size(gtBoxes{cI}, 1),1);
end

% indices per class
classIdx = ones(1, nClasses);

for cI = 1:length(gtBoxes)
    for i = 1:size(gtBoxes{cI}, 1)
        boScores{cI}(classIdx(cI)) = ...
            BoxBestOverlap(gtBoxes{cI}(i,:), testBoxes{gtNrs{cI}(i)});
        classIdx(cI) = classIdx(cI) + 1;
    end
end

% Calculation abo and mabo measures
abo = zeros(nClasses, 1);
for cI = 1:nClasses
    abo(cI) = mean(boScores{cI});
end
mabo = mean(abo);

% Calculation avgNumBoxes
numBoxes = zeros(length(testBoxes), 1);
for i=1:length(testBoxes)
    numBoxes(i) = size(testBoxes{i}, 1);
end
avgNumBoxes = mean(numBoxes);