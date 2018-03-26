function [ cimages ] = getSurfaceLabelInitial( img )
%UNTITLED given image, output surface label
%   Detailed explanation goes here
global config

bufdir = './Buffer/';
codedir = './Toolbox/SpatialLayout_shrink/';

% imdir = './roomModel/SpatialLayout_shrink/Images_resized/';
% moveFolder = './roomModel/SpatialLayout_shrink/';
% imsegdir = './roomModel/SpatialLayout_shrink/Imsegs/';

t= clock();
imagename = sprintf('%02d_%02d_%02d_%02d', t(4), t(5), int32(t(6)*10000), randi(99,1));
segext = '.pnm';
imgext = '.ppm';
imwrite( img, [bufdir imagename '.jpg']);
imwrite( img, [bufdir imagename imgext]);

bufImgName = [bufdir imagename imgext];
bufResName = [bufdir imagename segext];

% cmdLine = sprintf('D:\\YINDA\\SceneParsing\\rectangleDetector\\segmentation\\segment 0.8 100 100 %s %s', bufImgName, bufResName);
% cmdLine = sprintf('./rectangleDetector/segmentation/segment 0.8 100 100 %s %s', bufImgName, bufResName);
cmdLine = sprintf([config.segLocation '0.8 100 100 %s %s'], bufImgName, bufResName);
system(cmdLine);   %0.8, 100, 100


% fn=[moveFolder '/Imsegs/' imagename(1:end-4) '.' segext];
imseg = processSuperpixelImage(bufResName);

tic
imdata = mcmcComputeImageData(im2double(img), imseg);% made changes here
toc

% load(fullfile([moveFolder '/LabelClassifier/'], 'Classifiers_gc.mat'));
load([codedir '/LabelClassifier/Classifiers_gc.mat']);

spfeatures = mcmcGetAllSuperpixelData(bufdir, imseg);
[efeatures, adjlist] = mcmcGetAllEdgeData(spfeatures, imseg(1));

nsegments=[5 15 25 35 40 60 80 100];
pE{1} = test_boosted_dt_mc(eclassifier, efeatures{1});
pE{1} = 1 ./ (1+exp(ecal(1)*pE{1}+ecal(2)));
smaps{1} = msCreateMultipleSegmentations(pE{1}, adjlist{1}, ...
    imseg(1).nseg, nsegments);

for k = 1:numel(nsegments)
    if max(smaps{1}(:, k))>0
        segfeatures{1, k} = mcmcGetSegmentFeatures(imseg, ...
            spfeatures{1}, imdata, smaps{1}(:, k), (1:max(smaps{1}(:, k))));
    end
end

%Get surface label confidences initial from GC
normalize = 1;
pg=zeros(imseg.nseg,7);%7 labels
pg = msTest(imseg, segfeatures, smaps, ...
    labelclassifier, segclassifier,normalize);

cimages = msPg2confidenceImages(imseg,pg);
% sepScene(sid).cimages = cimages;
%     viewSurfaceLabel(img, cimages);
delete([bufdir imagename '.jpg']);
delete(bufImgName);
delete(bufResName);

end

