function [ labelMap ] = gbImageSegment( img, sigma, K, min )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

c = clock;
bufImgName = sprintf('./rectangleDetector/segmentation/img_%d_%d_%d.ppm', ...
                    c(end-2), c(end-1), c(end)*1000);
imwrite(img, bufImgName);
bufResName = sprintf('./rectangleDetector/segmentation/res_%d_%d_%d.ppm', ...
                    c(end-2), c(end-1), c(end)*1000);

cmdLine = sprintf('D:\\SceneParsing\\rectangleDetector\\segmentation\\segment %f %d %d %s %s', ...
                    sigma, K, min, bufImgName, bufResName);
system(cmdLine);

B = imread(bufResName);
C = unique(reshape(B,[],3),'rows');
labelMap = rgb2ind(B, double(C)./255);

delete( bufImgName);
delete( bufResName);

end

