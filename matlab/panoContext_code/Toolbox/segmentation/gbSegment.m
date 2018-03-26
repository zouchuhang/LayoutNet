function [ temp ] = gbSegment( img, sigma, k, minSz )
%GBSEGMENT Summary of this function goes here
%   Detailed explanation goes here
[height, width, ~] = size(img);
img_smooth = smooth(img, sigma);
edges = zeros(width*height*4,3);
[gridX, gridY] = meshgrid(1:width,1:height);
num = 0;

x = gridX(:);
y = gridY(:);
vector = [1 0; 0 1; 1 1; 1 -1];
inda = sub2ind([height width], y, x);

for vid = 1:size(vector,1)
    xv = x+vector(vid,1);
    yv = y+vector(vid,2);
    valid = xv>=1 & xv<=width & yv>=1 & yv<=height;
    indav = inda(valid);
    indbv = sub2ind([height width], yv(valid), xv(valid));
    diff = (img_smooth(indav)-img_smooth(indbv)).^2 ...
         + (img_smooth(indav+width*height)-img_smooth(indbv+width*height)).^2 ...
         + (img_smooth(indav+2*width*height)-img_smooth(indbv+2*width*height)).^2;
    edges(num+1:num+sum(valid),:) = [x(valid)-1 + (y(valid)-1)*width ...
                                     xv(valid)-1 + (yv(valid)-1)*width ...
                                     sqrt(diff)];
    num = num + sum(valid);
end

edges = edges';
segment = segmentGraphMex(width, height, num, edges, k, minSz);
L = unique(segment);

temp = zeros(size(segment));
for i = 1:length(L)
    temp(segment==L(i)) = i;
end


end

