function [ panoReport ] = viewRectangles( img, rect, num, scr )
%VIEWRECTANGLES Summary of this function goes here
%   Detailed explanation goes here
if isempty(num)
    num = sum(rect.score>scr);
else
    num = min(num, length(rect.score));
end

rects = zeros(0,12);
for j = 1:num
    rects = [rects; rect.xyzBox(j, :)];
end

lines = zeros(0,6);
for i = 1:size(rects,1)
    rect = [rects(i,1:3); rects(i,4:6); rects(i,7:9); rects(i,10:12)];
    lines = [lines; lineFromTwoPoint(rect([1 2 3 4],:), rect([2 3 4 1],:))];
end
panoReport = paintParameterLine( lines, 1024, 512, img);

end

