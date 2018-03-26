% mcmcTrainScript
% 5-fold cross-validation training

LOAD = 0;

ncv = 5;

imdir = '/IUS/vmr7/dhoiem/context/images/all_images';
outdir = '/IUS/vmr7/dhoiem/context/data';

if LOAD
    load([outdir '/rand_indices.mat']);
    load([outdir '/allimsegs2.mat']);
    load([outdir '/mcmcSuperpixelData.mat']);
    load([outdir '/mcmcSuperpixelClassifier.mat']);
    load([outdir '/mcmcEdgeData.mat']);
    load([outdir '/mcmcEdgeClassifier.mat']);
    load([outdir '/mcmcSegmentationClassifierCV.mat']);    
    load([outdir '/mcmcGtSegmentData2.mat']);
    load([outdir '/mcmcSegmentClassifierCV2.mat']);
end

nimages = numel(imsegs);

if ~exist('spfeatures')
    spfeatures = mcmcGetAllSuperpixelData(imdir, imsegs);
    save([outdir '/mcmcSuperpixelData.mat'], 'spfeatures');
end

if ~exist('vclassifierSP')
    [vclassifierSP, hclassifierSP] = ...
        mcmcTrainSuperpixelClassifier(spfeatures(cluster_images), imsegs(cluster_images));
    save([outdir '/mcmcSuperpixelClassifier.mat'], 'vclassifierSP', 'hclassifierSP');
end

if ~exist('efeatures')
    [efeatures, adjlist] = mcmcGetAllEdgeData(spfeatures, imsegs);
    save([outdir '/mcmcEdgeData.mat'], 'efeatures', 'adjlist');
end

if ~exist('eclassifier')
    eclassifier = mcmcTrainEdgeClassifier(efeatures(cluster_images), ...
        adjlist(cluster_images), imsegs(cluster_images));
    save([outdir '/mcmcEdgeClassifier.mat'], 'eclassifier');
end

% if ~exist('gtmaps')
%     [gtmaps, gtlabels] = mcmcGetGroundTruthSegmentations(imsegs, adjlist);
%     gtsegfeatures = mcmcGetAllSegmentFeatures(imsegs, imdir, gtmaps, spfeatures, vclassifierSP, hclassifierSP);
%     save([outdir '/mcmcGtSegmentPerfectData.mat'], 'gtmaps', 'gtlabels', 'gtsegfeatures');    
% end
% 
% if ~exist('vclassifier')
%     for k = 1:ncv       
%         trainind = setdiff([1:numel(imsegs)], ...
%             cv_images(floor((k-1)*numel(cv_images)/ncv)+1:floor(k*numel(cv_images)/ncv)));
%         [vclassifier(k), hclassifier(k)] = ...
%             mcmcTrainSegmentClassifier(imsegs(trainind), gtsegfeatures(trainind), gtlabels(trainind));
%     end
%     save([outdir '/mcmcSegmentClassifierPerfectCV.mat'], 'vclassifier', 'hclassifier');    
% end

if ~exist('segclassifier')
    [segdata, seglabels, segimind] = mcmcGenerateRandomSegments3(imsegs(cv_images), imdir,  ...
        adjlist(cv_images), spfeatures(cv_images), efeatures(cv_images), vclassifierSP, hclassifierSP, ...
        eclassifier);      
    
    for k = 1:ncv       
        trainind = setdiff([1:numel(cv_images)], ...
            [floor((k-1)*numel(cv_images)/ncv)+1:floor(k*numel(cv_images)/ncv)]);             
        segclassifier(k) = mcmcTrainSegmentationClassifier(segdata, seglabels, segimind, trainind);
    end
    save([outdir '/mcmcSegmentationClassifierCV.mat'], 'segclassifier', 'segdata', 'seglabels', 'segimind');
end

if ~exist('gtfeatures')
    [gtfeatures, gtlabels, gtlabelprobs, gtweights] = mcmcGetAllGoodSegments(...
        imsegs(cv_images), imdir, vclassifierSP, hclassifierSP, eclassifier, ...
        segclassifier, spfeatures(cv_images), adjlist(cv_images), efeatures(cv_images), ncv);   
    save([outdir '/mcmcGtSegmentData2.mat'], 'gtfeatures', 'gtlabels', 'gtlabelprobs', 'gtweights');
end

if ~exist('vclassifier')
    for k = 1:ncv       
        trainind = setdiff([1:numel(cv_images)], ...
            cv_images(floor((k-1)*numel(cv_images)/ncv)+1:floor(k*numel(cv_images)/ncv)));
        [vclassifier(k), hclassifier(k)] = ...
            mcmcTrainSegmentClassifier2(gtfeatures(trainind), gtlabels(trainind), gtweights(trainind));
    end
    save([outdir '/mcmcSegmentClassifierCV2.mat'], 'vclassifier', 'hclassifier');    
end

% if ~exist('vregressor')
%     for k = 1:ncv       
%         trainind = setdiff([1:numel(cv_images)], ...
%             cv_images(floor((k-1)*numel(cv_images)/ncv)+1:floor(k*numel(cv_images)/ncv)));
%         [vregressor(k), hregressor(k)] = ...
%             mcmcTrainSegmentRegressor(gtfeatures(trainind), gtlabelprobs(trainind), gtweights(trainind));
%     end
%     save([outdir '/mcmcSegmentRegressorCV.mat'], 'vregressor', 'hregressor');    
% end
