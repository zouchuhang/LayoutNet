function [ tempimg ] = visualCimage( img, cimages )
%VISUALCIMAGE Summary of this function goes here
%   Detailed explanation goes here

rl=[255  255   0  0    0    255  255];
gl=[255  0   255  255  0    0    255];
bl=[0    0   255  0    255  255  255];

[aa indd]=max(cimages,[],3);

% figure(102);
clear mask_color;
mask_r = rl(indd);
mask_g = gl(indd);
mask_b = bl(indd);
mask_color(:,:,1) = mask_r;
mask_color(:,:,2) = mask_g;
mask_color(:,:,3) = mask_b;

hsvmask=rgb2hsv(mask_color);
hsvmask(:,:,3)=aa*255;
%     hsvmask(:,:,2)=aa;
mask_color=hsv2rgb(hsvmask);

tempimg = double(img)*0.5 + mask_color*0.5;
%         tempimg =  mask_color;
% imshow(uint8(tempimg));

end

