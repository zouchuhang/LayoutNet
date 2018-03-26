function [ rectangles ] = pruneRectangle( rotImg, rectangles )
%PRUNERECTANGLE Prune rectangles based on gradients
%   Rectangles should present stronger gradients on edges

rotImg = imresize(rotImg, 0.25);
[H,W,~] = size(rotImg);
% img = zeros(H,W);

[FX, FY] = gradient(double(rgb2gray(rotImg)));
F = sqrt(FX.^2 + FY.^2);
[FD, ~, ~] = dt2(F, 0.1, 0, 0.1, 0 );

for vid = 1:6
    numrect = size(rectangles(vid).highPixelBox,1);
    valid = false(numrect, 1);
    rects = rectangles(vid).xyzBox;
    for i = 1:numrect
        rect = [rects(i,1:3); rects(i,4:6); rects(i,7:9); rects(i,10:12)];
        lines = lineFromTwoPoint(rect([1 2 3 4],:), rect([2 3 4 1],:));
        panoReport = paintParameterLine( lines, W, H);
        G = FD(panoReport(:)>0);
        SG = sort(G, 'descend');
        N = round(length(G)/2);
        
%         valid(i) = sum(SG(1:N))/N;
        if sum(SG(1:N))/N > 20
            valid(i) = true;
        end
    end
    
    rectangles(vid).highPixelBox = rectangles(vid).highPixelBox(valid,:);
    rectangles(vid).count = sum(valid);
    rectangles(vid).xyzBox = rectangles(vid).xyzBox(valid,:);
    rectangles(vid).score = rectangles(vid).score(valid);
end

end

