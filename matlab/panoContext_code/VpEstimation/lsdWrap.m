function [ edgeMap, edgeList ] = lsdWrap( img, cmdline )
%LSDWRAPPER Summary of this function goes here
%   Detailed explanation goes here
global config

if ~exist('cmdline', 'var')
    cmdline = '';
end
% rng('shuffle');
t= clock();
randomID = sprintf('%02d_%02d_%02d_%02d', randi(99,1), t(4), t(5), int32(t(6)*10000));
imgbuf = [config.bufferFolder 'buf_' randomID '.pgm'];
edgbuf = [config.bufferFolder 'buf_result_' randomID '.txt'];
%% line detection
while exist('imgbuf','file') || exist('edgbuf','file')
    fprintf('buf name conflict\n');
    pause(0.1);
end    
imwrite( img, imgbuf, 'pgm');
% system(['D:\YINDA\SceneParsing\edgeDetection\lsd ' cmdline ' ' imgbuf ' ' edgbuf]);
% system(['./edgeDetection/lsd_1.6/lsd ' cmdline ' ' imgbuf ' ' edgbuf]);
system([config.lsdLocation cmdline ' ' imgbuf ' ' edgbuf]);

edgeList = load( edgbuf);
edgeNum = size(edgeList,1);

%% draw on image
[imgH, imgW, ~] = size(img);
edgeMap = zeros(imgH, imgW);
for i = 1:edgeNum
    x = linspace( edgeList(i,1)+1, edgeList(i,3)+1, 1000);
    y = linspace( edgeList(i,2)+1, edgeList(i,4)+1, 1000);
    xx = max( min( round(x), imgW), 1);
    yy = max( min( round(y), imgH), 1);
    index = sub2ind( [imgH imgW], yy, xx);
    edgeMap(index) = 1;%rand(1);
end

%% delete buffer files
delete( imgbuf );
delete( edgbuf );
