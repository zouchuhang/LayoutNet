function [ panoLineImg ] = viewRegionHypothesis( img, AllCuboid, thres )
%VIEWREGIONHYPOTHESIS Summary of this function goes here
%   Detailed explanation goes here

AllLines = zeros(0,6);
for cid = 1:AllCuboid.count
    if AllCuboid.score(cid)<thres
        continue;
    end
        
    view = AllCuboid.views(cid,:);
    xyzBox = AllCuboid.xyzBox(cid,:);
    numRectangle = sum(view~=0);
    
    for rid = 1:numRectangle
        rect = reshape(xyzBox((rid*12-11):(rid*12)), 3, 4)';
        lines = lineFromTwoPoint(rect([1 2 3 4], :), rect([2 3 4 1], :));
        AllLines = [AllLines; lines];
    end
end

% img = imresize( im2double(img), [512 1024]);
% panoLineImg = paintParameterLine( AllLines, 1024, 512, img);
img = imresize( im2double(img), [1024 2048]);
panoLineImg = paintParameterLine( AllLines, 2048, 1024, img);

% imshow(panoLineImg);



end

