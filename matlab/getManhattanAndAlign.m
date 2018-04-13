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
[ olines, vp, views, edges, panoEdge, score, angle] = panoEdgeDetection( im_o, imgSize, qError);

% align image
[rotImg, R] = rotatePanorama(im_o, vp(3:-1:1,:));


