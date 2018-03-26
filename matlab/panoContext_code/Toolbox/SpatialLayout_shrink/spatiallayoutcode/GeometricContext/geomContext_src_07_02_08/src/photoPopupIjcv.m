function photoPopupIjcv(fnData, fnImage, extSuperpixelImage, outdir)
% photoPopup [fnData] [fnImage] [extSuperpixelImage] [outdir]
% fnData: filename for .mat file containing classifier data
% fnImage: filename for original RGB image
% extSuperpixelImage: filename extension for superpixel image 
% outdir: directory to output results
%
% Example usage: 
% Have data file classifiers_08_22_2005 stored in current directory, the
% original image ../images/103_0354.jpg, and the superpixel image
% '../images/103_0354.pnm.  I want results to go in ../results.  
%
% photoPopup ./classifiers_08_22_2005 ../images/103_0354.jpg pnm ../results 
% 
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

if nargin~=4
    disp(' ')
    disp(['Only ' num2str(nargin) ' arguments given; Requires 4'])
    disp(' ')
    disp('photoPopup [fnData] [fnImage] [extSuperpixelImage] [outdir]')
    disp('fnData: filename for .mat file containing classifier data')
    disp('fnImage: filename for original RGB image')
    disp('extSuperpixelImage: filename extension for superpixel image')
    disp('outdir: directory to output results')
    disp(' ')
    disp('Example usage:')
    disp('Have data file classifiers_08_22_2005 stored in current directory, the')
    disp(' original image ../images/103_0354.jpg, and the superpixel image');
    disp(' ../images/103_0354.pnm.  I want results to go in ../results. '); 
    disp(' ');
    disp('photoPopup ./classifiers_08_22_2005 ../images/103_0354.jpg pnm ../results ')
    disp(' ')
    disp('Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005')
    disp('Permission granted to non-commercial enterprises for ')
    disp('modification/redistribution under GNU GPL.')
    disp('Current Version: 1.0  09/30/2005')
    disp(' ')
    return;
end

% load classifiers
classifiers = load(fnData);

% read image
im = im2double(imread(fnImage));

% separate directory from filename
[imdir, fn] = strtok(fnImage, '/');
if isempty(fn)
    fn = fnImage;
    imdir = '.';
else
    while ~strcmp(strtok(fn, '/'), fn(2:end))
        [rem, fn] = strtok(fn, '/');
        imdir = [imdir '/' rem];
    end
    fn = fn(2:end);
end

if isempty(extSuperpixelImage)
    
    imsegs = im2superpixels(im);
    imsegs.imname = fn;
else

    % compute superpixel map
    imsegs.image_name = fn;
    imsegs = APPimages2superpixels(imdir, extSuperpixelImage, imsegs);

end
    
% compute geometry labels and confidences            
pg = ijcvTestImage(im, imsegs, classifiers);               

% get confidence maps
[bn, ext] = strtok(fn, '.');
% [cimages, cnames] = APPclassifierOutput2confidenceImages(imsegs, conf_map);
% for i = 1:length(cnames)
%     imwrite(cimages{1}(:, :, i), [outdir '/' bn '.' cnames{i} '.pgm']);
% end
% 
% % get labeled image
% limage = APPgetLabeledImage(im, imsegs, labels.vert_labels, labels.vert_conf, ...
%         labels.horz_labels, labels.horz_conf);     
% imwrite(limage, [outdir, '/', bn, '.l.jpg']);
% imwrite(im, [outdir, '/', fn]);
[pv, ph] = splitpg(pg);
lim = APPgetLabeledImage2(im, imsegs, pv, ph);
imwrite(lim, [outdir, '/', bn, '.l.jpg']);
imwrite(im, [outdir, '/', fn]);

% get popup model
lines = APPgetLargeConnectedEdges(rgb2gray(im), min([size(im, 1) size(im, 2)]*0.02), imsegs);
hy = 1-APPestimateHorizon(lines);
[maxval, vlabels] = max(pv, [], 2);
APPwriteVrmlModel_v2(imdir, imsegs, vlabels, hy, outdir)   
   
