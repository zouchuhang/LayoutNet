function updateSSHypothesis( aid )
%UPDATESSHYPOTHESIS Summary of this function goes here
%   Detailed explanation goes here
%matlabpool(2);
[~,b] = system('hostname');
fprintf('Running on %s\n',b);

% diaryfilename = sprintf('../diary/diary%03d.txt', aid);
% diary(diaryfilename);
% diary on;

addpath(genpath('/n/fs/sun360/panoParser/SceneParsing'));
cd '/n/fs/sun360/panoParser/SceneParsing';
folderName = './annotation/dataset/bedroom/';
% folderName = './annotation/dataset/living_room/';

load([folderName '/IMGLIST.mat']);
load([folderName '/ANNO_ALL.mat']);
fprintf('Data loaded\n');
testfile2 = [folderName IMGLIST(aid).name '/bottomupcuboid_P1.mat'];
if exist(testfile2, 'file')
	fprintf('Job: %d: Previously done^_^\n', aid);
	return;
end

anno_valid = false( 1, length(ANNO_ALL));
for i = 1:length(anno_valid)
    anno_valid(i) = ~isempty(ANNO_ALL(i).ANNO3D) && ANNO_ALL(i).ANNO3D.b_singleroom;
end

if ~anno_valid(aid)
    return;
end

rotImg = imread([IMGLIST(aid).small_imgname(1:end-9) 'rotate.png']);
fprintf('Image loaded\n');

SS = getSelectiveSearch(rotImg);

if exist(testfile2, 'file')
	fprintf('Job: %d: Previously done^_^\n', aid);
	return;
end

AllCuboid = regionBasedHypothesisFromSS(SS, 1024, 2048);

if exist(testfile2, 'file')
	fprintf('Job: %d: Previously done^_^\n', aid);
	return;
end

rectangle = rectDetectionFlexFromSS(SS, 1024, 2048);

if exist(testfile2, 'file')
	fprintf('Job: %d: Previously done^_^\n', aid);
	return;
end

AllCuboid_NEW = pruneCuboid(rotImg, AllCuboid);
rectangle_NEW = pruneRectangle(rotImg, rectangle);

load([folderName IMGLIST(aid).name '/bottomupcuboid.mat']);
fprintf('Old data loaded\n');

regionCuboid.views = [regionCuboid.views; AllCuboid_NEW.views];
regionCuboid.xyzBox = [regionCuboid.xyzBox; AllCuboid_NEW.xyzBox];
regionCuboid.score = [regionCuboid.score; AllCuboid_NEW.score];
regionCuboid.count = regionCuboid.count + AllCuboid_NEW.count;

for i = 1:6
    rectangle(i).highPixelBox = [rectangle(i).highPixelBox; rectangle_NEW(i).highPixelBox];
%     rectangle(i).count = rectangle(i).count + rectangle_NEW(i).count;
    rectangle(i).xyzBox = [rectangle(i).xyzBox; rectangle_NEW(i).xyzBox];
    rectangle(i).score = [rectangle(i).score; rectangle_NEW(i).score];
end

save([folderName IMGLIST(aid).name '/bottomupcuboid_P1.mat'], ...
            'rectangle', 'rectCuboid2', 'rectCuboid3', 'regionCuboid');
fprintf('save data for bottomcuboid_P1.mat\n');

%% move updatePruneRectDetBySegmentation
testfile2 = [folderName IMGLIST(aid).name '/rectangleSegScore_P1.mat'];
[ rect_score ] = getSegScrOfRect( rotImg, rectangle );
save(testfile2, 'rect_score');
fprintf('save data for rectangleSegScore_P1.mat\n');

end

