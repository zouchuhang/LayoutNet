% trainScript3: NO LONGER IN USE
 
imdir = '../images/all_images';
outdir = '/IUS/vmr20/dhoiem/data/mcmcdata';
ncv = 5;

nsegments = [3 5 6 6 7 8 9 10 11 12 13 14 15 17 18 20 23 27 32 45 50 75 100];
%nsegments = [3 4 5 6 7 7 8 8 9 9 10 11 11 12 12 13 13 14 14 15 16 16 17 19 20 21 23 25 28 31 37 51];

if ~exist('labdata')
    % gather data
    for tf = 1:numel(cv_images)

        f = cv_images(tf);

        disp([num2str(tf) ': ' imsegs(f).imname])

        [pvSP{tf}, phSP{tf}, pE{tf}] = mcmcInitialize(spfeatures{f}, efeatures{f}, ...
            adjlist{f}, imsegs(f), vclassifierSP, hclassifierSP, eclassifier, 'labels');
        smaps{tf} = generateMultipleSegmentations(pE{tf}, adjlist{f}, imsegs(f).nseg, nsegments);

        im = im2double(imread([imdir '/' imsegs(f).imname]));
        imdata = mcmcComputeImageData(im, imsegs(f));
%         imdata.pvSP = pvSP;
%         imdata.phSP = phSP;   

        for k = 1:numel(nsegments)
            labdata{tf, k} = mcmcGetSegmentFeatures(imsegs(f), spfeatures{f}, imdata, smaps{tf}(:, k), (1:max(smaps{tf}(:, k))));
            segdata{tf, k} = mcmcGetSegmentationFeatures(pvSP{tf}, phSP{tf}, pE{tf}, adjlist{f}, imsegs(f).npixels, smaps{tf}(:, k), (1:max(smaps{tf}(:, k))));
        end
    end
    [vlab, hlab, vw, hw] = segmentation2labels2(imsegs(cv_images), smaps);
    save([outdir '/traindata5.mat'], 'smaps', 'labdata', 'segdata', 'mclab', 'mcprc', 'allprc', 'seglabel', 'unilabel', 'trainw', 'pvSP', 'phSP', 'pE');
end

for k = 1:ncv
    disp(['Iteration: ' num2str(k)]);
    testind{k} = (floor((k-1)*numel(cv_images)/ncv)+1):(floor(k*numel(cv_images)/ncv));
    trainind{k} = setdiff([1:numel(cv_images)], testind{k});
    svclassifier(k) = mcmcTrainSegmentationClassifier2(segdata(trainind{k}, :), vlab(trainind{k}, :), vw(trainind{k}, :));    
    shclassifier(k) = mcmcTrainSegmentationClassifier2(segdata(trainind{k}, :), hlab(trainind{k}, :), hw(trainind{k}, :));    
    [vclassifier(k), hclassifier(k)] = ...
        mcmcTrainSegmentClassifier3(labdata(trainind{k}, :), vlab(trainind{k}, :), hlab(trainind{k}, :), vw(trainind{k}, :), hw(trainind{k}, :)); 
end
[vacc, hacc, vcm, hcm] = testMultipleSegmentationsCV3(imsegs(cv_images), ...
    labdata, segdata, smaps, vclassifier, hclassifier, svclassifier, shclassifier, pvSP, phSP, ncv);

save([outdir '/results5b.mat'], 'vclassifier', 'hclassifier', 'svclassifier', 'shclassifier', 'vacc', 'hacc', 'vcm', 'hcm');
