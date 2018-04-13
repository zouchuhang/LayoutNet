% path
addpath(genpath('./panoContext_code/'));

% param
imgSize = 320;
qError = 0.7;

% read image
im = imread('pano_arrsorvpjptpii.jpg');

global config;
config.NewProjFolder = './panoContext_code/';
config.lsdLocation = [config.NewProjFolder 'VpEstimation/lsd '];
config.bufferFolder = './Buffer_gt/';
if ~exist(config.bufferFolder, 'dir')
    mkdir(config.bufferFolder);
end

[im_h, im_w, im_dim] = size(im);
im_o = im2double(im);

% extract manhattan line
[ olines, vp, views, edges, panoEdge, score, angle] = panoEdgeDetection_line( im_o, imgSize, qError);
panoEdge = double(panoEdge(:,:,[1,4,7])>0);

% align image
[rotEdge, R] = rotatePanorama(panoEdge, vp(3:-1:1,:));
Img_small = imresize(im_o, [1024 2048]);
[rotImg, R] = rotatePanorama(Img_small, vp(3:-1:1,:));
panoImg = imresize(rotImg, [im_h im_w]);

% visual
subplot(2,1,1)
imagesc(rotImg); axis image
subplot(2,1,2)
imagesc(rotEdge); axis image

