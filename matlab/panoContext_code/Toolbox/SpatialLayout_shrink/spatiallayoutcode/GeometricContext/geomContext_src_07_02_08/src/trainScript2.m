% trainScript2
 
imdir = '../images/all_images';
outdir = '/IUS/vmr20/dhoiem/data/mcmcdata';
ncv = 5;

nsegments = [5:5:50 60:10:100];
%nsegments = [3 5 7 9 11 13 15 20 25 30 35 40 50 60 70 80 90 100];
%nsegments = [3 5 6 6 7 8 9 10 11 12 13 14 15 17 18 20 23 27 32 45 50 75 100];
%nsegments = [3 4 5 6 7 7 8 8 9 9 10 11 11 12 12 13 13 14 14 15 16 16 17 19 20 21 23 25 28 31 37 51];

% Computes the superpixel features
if ~exist('spfeatures')
    spfeatures = mcmcGetAllSuperpixelData(imdir, imsegs);
    save([datadir '/superpixelData.mat'], 'spfeatures');
end

% Trains the superpixel classifiers
if ~exist('vclassifierSP')
    trainind = cluster_images;
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
    trainind = cluster_images;
    eclassifier = mcmcTrainEdgeClassifier(efeatures(trainind), ...
        adjlist(trainind), imsegs(trainind));
    ecal = calibrateEdgeClassifier(efeatures(trainind), adjlist(trainind), ...
        imsegs(trainind), eclassifier, 1);
    ecal = ecal{1};
    save([datadir '/edgeClassifier.mat'], 'eclassifier', 'ecal');
end


if ~exist('labdata')
    % gather data
    for tf = 1:numel(cv_images)

        f = cv_images(tf);                
        
        disp([num2str(tf) ': ' imsegs(f).imname])

        [pvSP{tf}, phSP{tf}, pE{tf}] = mcmcInitialize(spfeatures{f}, efeatures{f}, ...
            adjlist{f}, imsegs(f), vclassifierSP, hclassifierSP, eclassifier, 'none');
        smaps{tf} = generateMultipleSegmentations2(pE{tf}, adjlist{f}, imsegs(f).nseg, nsegments);

        im = im2double(imread([imdir '/' imsegs(f).imname]));
        imdata = mcmcComputeImageData(im, imsegs(f));
%         imdata.pvSP = pvSP;
%         imdata.phSP = phSP;   

        for k = 1:numel(nsegments)
            labdata{tf, k} = mcmcGetSegmentFeatures(imsegs(f), spfeatures{f}, imdata, smaps{tf}(:, k), (1:max(smaps{tf}(:, k))));
            segdata{tf, k} = mcmcGetSegmentationFeatures(pvSP{tf}, phSP{tf}, pE{tf}, adjlist{f}, imsegs(f).npixels, smaps{tf}(:, k), (1:max(smaps{tf}(:, k))));
            [mclab{tf, k}, mcprc{tf, k}, allprc{tf, k}, trainw{tf, k}] = segmentation2labels(imsegs(f), smaps{tf}(:, k));
            unilabel{tf, k} = mclab{tf, k}.*(mcprc{tf, k}>0.99);
            seglabel{tf,k} =  1*(mcprc{tf, k}>0.99) + (-1)*(mcprc{tf, k}<0.95);
            
            if k==8
                smapk = smaps{tf}(:, k);
                figure(1), imagesc(label2rgb(smapk(imsegs(f).segimage))), axis image
                drawnow; pause(1);
            end            
            
        end
    end
    save([outdir '/traindata7.mat'], 'smaps', 'labdata', 'segdata', 'mclab', 'mcprc', 'allprc', 'seglabel', 'unilabel', 'trainw', 'pvSP', 'phSP', 'pE');
end

for k = 1:ncv
    disp(['Iteration: ' num2str(k)]);
    testind{k} = (floor((k-1)*numel(cv_images)/ncv)+1):(floor(k*numel(cv_images)/ncv));
    trainind{k} = setdiff([1:numel(cv_images)], testind{k});
    sclassifier(k) = mcmcTrainSegmentationClassifier2(segdata(trainind{k}, :), seglabel(trainind{k}, :), trainw(trainind{k}, :));    
    [vclassifier(k), hclassifier(k)] = ...
        mcmcTrainSegmentClassifier2(labdata(trainind{k}, :), unilabel(trainind{k}, :), trainw(trainind{k}, :));
    %[vregressor(k), hregressor(k)] = ...
    %    trainRegressorLR(labdata(trainind, :), allprc(trainind, :), trainw(trainind, :));     
    %[sregressor(k)] = ...
    %    trainRegressorLR(labdata(trainind{k}, :), mcprc(trainind{k}, :), trainw(trainind{k}, :));         
end
[vacc, hacc, vcm, hcm, pg] = testMultipleSegmentationsCV2(imsegs(cv_images), ...
    labdata, segdata, smaps, vclassifier, hclassifier, sclassifier, pvSP, phSP, ncv);

save([outdir '/results7.mat'], 'vclassifier', 'hclassifier', 'sclassifier', 'vacc', 'hacc', 'vcm', 'hcm', 'pg');
