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
    load(fullfile(datadir, 'msClassifiers_tu_tmp.mat'));
    load(fullfile(datadir, 'msSegmentFeatures_tu_tmp.mat'));
    load(fullfile(outdir, 'msMultipleSegmentations_tu_tmp.mat'));    
end


nimages = numel(imsegs);


% disp('Converting labels to superpixels')
% if ~exist('splabels_test', 'var')
%     splabels_test = msLabelMap2Sp(labels, {imsegs.segimage});
%     save(fullfile(outdir, 'msSuperpixelLabelsTest.mat'), 'splabels_test');
% end

disp('Getting accuracy and confusion matrix')
if ~exist('acc', 'var')
    normalize = 0;
    pg = msTest(imsegs(test), segfeatures(test, :), smaps(test), ...
        labelclassifier, segclassifier, normalize);
    ignore = [5 8];
    [acc, cm, classcount] = msAnalyzeResult(imsegs(test), labels(test), pg, 0, ignore);
    cm(ignore, :) = []; cm(:, ignore) = [];        
    save(fullfile(outdir, 'msResult_tu_tmp.mat'), 'acc', 'cm', 'pg', 'testind', 'classcount', 'ignore');
    disp(['Accuracy: ' num2str(acc)])
    showConfusionMatrix(cm, classnames(setdiff(2:24, ignore+1)), 1)
end
