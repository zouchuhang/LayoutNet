% trainTestSurfaceContextScript
% 
% Assumes you have created:
%   imsegs  :  data structure that contains information on all training
%              and testing images
%   trainind:  the indices of the training images
%   testind:   the indices of the testing images

basedir = '???'; % set this appropriately
imdir = [basedir '/images/all_images'];
datadir = [basedir '/data'];
resultsdir = [basedir '/results'];

DO_TRAIN = 1;
DO_TEST = 1;

if DO_TRAIN

    nsegments = [10 20 30 40 50 60 80 100];

    % Computes the superpixel features
    if ~exist('spfeatures')
        spfeatures = mcmcGetAllSuperpixelData(imdir, imsegs);
        save([datadir '/superpixelData.mat'], 'spfeatures');
    end

    % Trains the superpixel classifiers
    if ~exist('vclassifierSP')
        [vclassifierSP, hclassifierSP] = ...
            mcmcTrainSuperpixelClassifier(spfeatures(trainind), imsegs(trainind));
        save([datadir '/superpixelClassifier.mat'], 'vclassifierSP', 'hclassifierSP');
    end

    % Computes the same-label features
    if ~exist('efeatures')
        [efeatures, adjlist] = mcmcGetAllEdgeData(spfeatures, imsegs);
        save([datadir '/edgeData.mat'], 'efeatures', 'adjlist');
    end

    % Trains the same-label classifier
    if ~exist('eclassifier')        
        eclassifier = mcmcTrainEdgeClassifier(efeatures(trainind), ...
            adjlist(trainind), imsegs(trainind));
        ecal = calibrateEdgeClassifier(efeatures(trainind), adjlist(trainind), ...
            imsegs(trainind), eclassifier, 1);
        ecal = ecal{1};
        save([datadir '/edgeClassifier.mat'], 'eclassifier', 'ecal');
    end
    
    % Computes the multiple segmentations, the segment features, and the
    % ground truth labels for each segment
    if ~exist('labdata')
        % gather data
        for f = 1:numel(imsegs)

            disp([num2str(f) ': ' imsegs(f).imname])

            [pvSP{f}, phSP{f}, pE{f}] = mcmcInitialize(spfeatures{f}, efeatures{f}, ...
                adjlist{f}, imsegs(f), vclassifierSP, hclassifierSP, eclassifier, ecal, 'none');
            smaps{f} = generateMultipleSegmentations2(pE{f}, adjlist{f}, imsegs(f).nseg, nsegments);

            im = im2double(imread([imdir '/' imsegs(f).imname]));
            imdata = mcmcComputeImageData(im, imsegs(f));

            for k = 1:numel(nsegments)
                labdata{f, k} = mcmcGetSegmentFeatures(imsegs(f), spfeatures{f}, imdata, smaps{f}(:, k), (1:max(smaps{f}(:, k))));
                [mclab{f, k}, mcprc{f, k}, allprc{f, k}, trainw{f, k}] = segmentation2labels(imsegs(f), smaps{f}(:, k));
                unilabel{f, k} = mclab{f, k}.*(mcprc{f, k}>0.95);
                seglabel{f,k} =  1*(mcprc{f, k}>0.95) + (-1)*(mcprc{f, k}<0.95);                                  
            end
        end
        save([datadir '/allData.mat'], 'smaps', 'labdata', 'segdata', 'mclab', 'mcprc', 'allprc', 'seglabel', 'unilabel', 'trainw', 'pvSP', 'phSP', 'pE');
    end
        
    % Trains the segment classifiers
    if ~exist('vclassifier')     
        sclassifier = mcmcTrainSegmentationClassifier2(labdata(trainind, :), seglabel(trainind{k}, :), trainw(trainind, :)); 
        [vclassifier, hclassifier] = ...
            mcmcTrainSegmentClassifier2(labdata(trainind, :), unilabel(trainind, :), trainw(trainind, :), 50000);             
    end

    save([datadir '/allClassifiers.mat'], 'vclassifier', 'hclassifier', 'sclassifier', 'eclassifier', 'vclassifierSP', 'hclassifierSP');

end

if DO_TEST
    
    % Computes the label confidences for each superpixel in the test images
    % and gives the final accuracy.
    % pg{image number}(superpixel number, [000 left center right porous solid sky])
    % gives the superpixel label confidences                                  
    [vacc, hacc, vcm, hcm, pg] = ...
        testMultipleSegmentationsCV2(imsegs(testind), labdata(testind, :), ...
        labdata(testind, :), smaps(testind), ...
        vclassifier, hclassifier, sclassifier, pvSP(testind), phSP(testind), 1);
    save([datadir '/results.mat'], 'vacc', 'hacc', 'vcm', 'hcm', 'pg');
    
    % Computes and writes the labeled images
    for f = 1:testind
        im = im2double(imread([imdir '/' imsegs(f).imname]));
        [pv, ph] = splitpg(pg{f});
        lim = APPgetLabeledImage2(im, imsegs(f), pv, ph);
        imwrite(im, [resultsdir '/' imsegs(f).imname]);
        imwrite(lim, [resultsdir '/' strtok(imsegs(f).imname, '.') '.l.jpg']);
    end
end