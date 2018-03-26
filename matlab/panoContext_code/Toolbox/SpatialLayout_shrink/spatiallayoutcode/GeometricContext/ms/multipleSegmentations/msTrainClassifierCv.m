% msTrainClassifiersCv

LOAD = 1;

nsegments = [5 7 10 15 20 25 35 50 60 70 80 90 100 Inf]; % number of segments per segmentation
labeltol = 0.9; % required percentage of single-label pixels for segment to be good
nclasses = 23;
ncv = 5;

datadir = '~/data/eccv08/msrc21/';
imdir = '~/data/msrc/MSRC_ObjCategImageDatabase_v2/Images/';
outdir = '~/data/eccv08/msrc21/';

if LOAD    
    load(fullfile(datadir, 'msrc_imsegs.mat'));
    load(fullfile(datadir, 'trainLabels.mat'));
    load(fullfile(datadir, 'trainTestFn_tu.mat'));
    load(fullfile(datadir, 'msSegmentFeatures_tu.mat'));
end

disp('Getting cross-validation subsets');
if ~exist('testcv', 'var')
    testcv = cell(ncv, 1);
    traincv = cell(ncv, 1);
    ind = train(randperm(numel(train)));
    fi = floor(numel(train)*(0:ncv-1)/ncv+1);
    li = floor(numel(train)*(1:ncv)/ncv);
    for k = 1:ncv
        testcv{k} = ind(fi(k):li(k));
        traincv{k} = setdiff(train, testcv{k});
    end
    save(fullfile(outdir, 'msCvind_tu.mat'), 'traincv', 'testcv');   
end
    
disp('Training segmentation classifier')
if ~exist('segclassifier', 'var')
    for k = 1:ncv
        disp(num2str(k))
        trainind = traincv{k};
        traindata = segfeatures(trainind, :);
        trainlabels = seggood(trainind, :);
        trainweights = trainw(trainind, :);
        segclassifier(k) = mcmcTrainSegmentationClassifier2(traindata, trainlabels, trainweights); 
    end
    save(fullfile(outdir, 'msSegmentationClassifier_cv_tu.mat'), 'segclassifier');    
end

disp('Training label classifier')
if ~exist('labelclassifier', 'var')
    for k = 1:ncv
        disp(num2str(k))
        trainind = traincv{k};
        traindata = segfeatures(trainind, :);
        trainlabels = seglabel(trainind, :);    
        trainweights = trainw(trainind, :);
        labelclassifier(k) = msTrainLabelClassifier(traindata, trainlabels, trainweights, ...
		classnames, Inf);
        %labelclassifier = mcmcTrainSegmentClassifier2(traindata, trainlabels, trainweights);   
    end
    save(fullfile(outdir, 'msLabelClassifier_cv_tu.mat'), 'labelclassifier');  
end

