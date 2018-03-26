function imsegs = APPimages2superpixels(input_dir, ext, gt)
% imsegs = APPimages2superpixels(input_dir, ext, gt)
% Create the imsegs structure from a segmentation image
%
% INPUT: 
% input_dir - The directory of segmentation images. Segments are denoted by
% different colors.
% ext - the extension of the segmentation image filenames
% gt - an existing imsegs structure (or [] if none exists)
% OUTPUT:
% imsegs - image segmentation data (and sometimes ground truth)
%
% Example: imsegs = APPimages2superpixels('../images', 'pnm', [])
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

if isempty(gt)
    files = dir([input_dir '/*.' ext]);
    for f = 1:length(files)
        gt(f).image_name = files(f).name;
    end
end
           

imsegs(length(gt)) = struct('imname', '', 'imsize', [0 0]);
for f = 1:length(gt)
    imname = gt(f).image_name;
    disp(imname)
    basename = strtok(imname, '.');
    im = imread([input_dir '/' basename '.' ext]);
    
    im = double(im);
    
    imsegs(f).imname = [strtok(imname,'.') '.jpg'];
    imsegs(f).imsize = size(im);
    imsegs(f).imsize = imsegs(f).imsize(1:2);
    im = im(:, :, 1) + im(:, :, 2)*256 + im(:, :, 3)*256^2;
    [gid, gn] = grp2idx(im);
    imsegs(f).segimage = uint16(reshape(gid, imsegs(f).imsize));
    imsegs(f).nseg = length(gn);

end
imsegs = APPgetSpStats(imsegs);
