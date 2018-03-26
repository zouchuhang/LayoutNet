function imsegs = processSuperpixelImage(fn)
% imsegs = processSuperpixelImage(fn)
% Creates the imsegs structure from a segmentation image
%
% INPUT: 
% fn - filenames of segmentation images. Use '/' (not '\') to separate directories. 
% Segments are denoted by different RGB colors.  
%
% OUTPUT:
% imsegs - image segmentation data 
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2006
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  04/24/2006
          
if isstr(fn)
    fn = {fn};
end

imsegs(length(fn)) = struct('imname', '', 'imsize', [0 0]);
for f = 1:length(fn)    
    toks = strtokAll(fn{f}, '/');
    imname = toks{end};
    basename = strtok(imname, '.');
    im = imread(fn{f});    
    im = double(im);
    
    imsegs(f).imname = [strtok(imname,'.') '.jpg'];
    imsegs(f).imsize = size(im);
    imsegs(f).imsize = imsegs(f).imsize(1:2);
    im = im(:, :, 1) + im(:, :, 2)*256 + im(:, :, 3)*256^2;
    [gid, gn] = grp2idx(im(:));
    imsegs(f).segimage = uint16(reshape(gid, imsegs(f).imsize));
    imsegs(f).nseg = length(gn);
end
imsegs = APPgetSpStats(imsegs);