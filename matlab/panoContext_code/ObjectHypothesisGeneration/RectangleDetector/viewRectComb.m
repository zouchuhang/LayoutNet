function [ panoReport ] = viewRectComb( img, cuboid, rectangle, r, num )
%VIEWRECTCOMB Summary of this function goes here
%   Detailed explanation goes here
num = min(num, size(cuboid(r).combo,1));
views = cuboid(r).views;
views(views==0) = [];

rects = zeros(0,12);
for j = 1:num
    combo = cuboid(r).combo(j,:);
    for i = 1:length(views)
        rects = [rects; rectangle(views(i)).xyzBox(combo(i), :)];
    end
end

lines = zeros(0,6);
for i = 1:size(rects,1)
    rect = [rects(i,1:3); rects(i,4:6); rects(i,7:9); rects(i,10:12)];
    lines = [lines; lineFromTwoPoint(rect([1 2 3 4],:), rect([2 3 4 1],:))];
end
panoReport = paintParameterLine( lines, 1024, 512, img);
% imshow(panoReport);

end

