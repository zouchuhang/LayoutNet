LOAD = 1;

ncv = 1; % number of cross-validation sets
nclasses = 23;

datadir = '~/data/eccv08/msrc21/';
imdir = '~/data/msrc/MSRC_ObjCategImageDatabase_v2/Images/';
outdir = '~/data/eccv08/msrc21/';

if LOAD    
    load(fullfile(datadir, 'msrc_imsegs.mat'));
    load(fullfile(datadir, 'testLabels.mat'));
    load(fullfile(datadir, 'trainTestFn_tu.mat'));    
    load(fullfile(datadir, 'msSegmentFeatures_tu.mat'));
    load(fullfile(outdir, 'msMultipleSegmentations_tu.mat'));
    load(fullfile(datadir, 'msSegmentationClassifier_cv_tu.mat'));
    load(fullfile(datadir, 'msLabelClassifier_cv_tu.mat'));
    load(fullfile(datadir, 'msCvind_tu.mat'));    
end

nimages = numel(imsegs);


% disp('Converting labels to superpixels')
% if ~exist('splabels_test', 'var')
%     splabels_test = msLabelMap2Sp(labels, {imsegs.segimage});
%     save(fullfile(outdir, 'msSuperpixelLabelsTest.mat'), 'splabels_test');
% end

disp('Getting accuracy and confusion matrix')
if ~exist('pg', 'var')
    
    for k = 1:numel(testcv)
        disp(num2str(k))
        classnums = [];
        for k2 = 1:numel(labelclassifier(k).names)
            classnums(k2) = str2num(labelclassifier(k).names{k2});
        end
        disp(['classnums: ' num2str(classnums)])
        
        ind = testcv{k};
        pg(ind) = msTest(imsegs(ind), segfeatures(ind, :), smaps(ind), ...
            labelclassifier(k), segclassifier(k), 0);                
        
        for k2 = 1:numel(ind)
            pg{ind(k2)}(:, classnums) = pg{ind(k2)};
        end
        
    end
%     ind = testcv{4};
%     for k = ind(:)'
%         pg{k} = [pg{k}(:, 1:4) zeros(size(pg{k}(:, 1))) pg{k}(:, 5:end)];
%     end

    pg = pg(train);
    
    
    ignore = [5 8];
    [acc, cm, classcount] = msAnalyzeResult(imsegs(train), labels(train), pg, 0, ignore);
    cm(ignore, :) = []; cm(:, ignore) = [];                

    showConfusionMatrix(cm, classnames(setdiff(2:24, ignore+1)), 1)    
    disp(['Accuracy: ' num2str(acc)])
    testind = train;
    save(fullfile(outdir, 'msResultCv_tu.mat'), 'acc', 'cm', 'pg', 'testind', 'classcount', 'ignore');
end
