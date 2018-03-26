function [labels, conf_map] = ...
    APPtestDirectory(segDensity, vClassifier, hClassifier, ...
        imdir, imsegs, varargin)
% [labels, conf_map] = APPtestDirectory(segDensity, vClassifier,
%                                   hClassifier, imdir, imsegs, varargin)
%
% Gets the geometry for each image with superpixels given by imsegs.
%
% Input:
%   segDensity: structure giving probability of 2 sp having same label
%   vClassifier: segment classifier for ground/vert/sky
%   hClassifier: segment classifier for subclassses of vert
%   imdir: the directory of the images referred to by imsegs
%   imsegs: the superpixel structure for each image
%   varargin{1} (optional): the output directory for displaying results
% Output:
%   labels: structure containing labeling results for each image
%   conf_map: the likelihoods for each sp for each class for each image
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

DO_PARALLEL = 0; % for running multiple parallel processes on directory

if length(varargin) > 0
    outdir = varargin{1};
end

for f = 1:length(imsegs)
           
    fn = imsegs(f).imname;
    bn = strtok(fn, '.');  
    
    if ~DO_PARALLEL || ~exist('outdir') || ~exist([outdir '/' bn '.c.mat'])
        
        if exist('outdir') && DO_PARALLEL % to mark as being processed
            save([outdir '/' bn '.c.mat'], 'fn');
        end        
            
        image = im2double(imread([imdir, '/', fn]));

        if size(image, 3) == 3 

            disp(['processing image ' fn]);

            [labels(f), conf_map(f), maps{f}, pmaps(f)] = ...
                    APPtestImage(image, imsegs(f), vClassifier, hClassifier, segDensity);               
            % for generating context
            [cimages, cnames] = APPclassifierOutput2confidenceImages(imsegs(f), conf_map(f));


            if length(varargin) > 0
                outdir = varargin{1};

                limage = APPgetLabeledImage(image, imsegs(f), labels(f).vert_labels, labels(f).vert_conf, ...
                    labels(f).horz_labels, labels(f).horz_conf);          
                imwrite(limage, [outdir, '/', bn, '.l.jpg']);
                imwrite(image, [outdir, '/', fn]);

                % for generating context
                glabels = labels(f);
                gconf_map = conf_map(f);
                gmaps = maps{f};
                gpmaps = pmaps(f);
                save([outdir '/' bn '.c.mat'], 'glabels', 'gconf_map', 'cimages', 'cnames', 'gmaps', 'gpmaps');

                % make a vrml file
                if 0
                    APPwriteVrmlModel(imdir, imsegs(f), labels(f), outdir);    
                end
            %else
            %    disp('warning: set up for ICCV, no output ==> no results')
            end

        end

        pause(0.05);

    end
    
end

for f = 1:length(imsegs)
    labels(f).imname = imsegs(f).imname;
end