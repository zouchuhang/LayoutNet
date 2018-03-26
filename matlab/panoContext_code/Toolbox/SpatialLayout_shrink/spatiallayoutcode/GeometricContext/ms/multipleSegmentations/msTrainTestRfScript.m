DO_TRAIN = 0;
DO_TEST = 1;



%% Training

if DO_TRAIN

datadir = '~/data/eccv08/msrc21/';

if 0     
    load(fullfile(datadir, 'msSegmentFeatures_tu_tmp.mat'));    
    load(fullfile(datadir, 'trainTestFn_tu.mat')); 
    %train = train(randperm(numel(train)));
end

maxdata = 1000000;
maxtree = 1000;
bagsize = 25000;

% good vs bad segment classifier
if 1
    
[data, lab, w] = msFormatClassifierData(...
    segfeatures(train, :), seggood(train, :), trainw(train, :), maxdata);
nfs = [1 2 4 8 16 32 48 64 94];
valset = (ceil(0.1*(numel(lab))):ceil(0.4*numel(lab)));
ntrees = 100;
bagfract = bagsize / (numel(lab)-numel(valset));
[segnf, segerr] = rfSelectPoolSize(data, (lab+1)/2+1, bagfract, ntrees, ...
    maxtree, nfs, valset, w);
end
ntrees = 200;
[dt_seg, segoob] = rfTrainForest(data, (lab+1)/2+1, ...
    bagfract, ntrees, maxtree, segnf, [], w);

% label classifier
[data, lab, w] = msFormatClassifierData(...
    segfeatures(train, :), seglabel(train, :), trainw(train, :), maxdata);
%valset = (1:ceil(0.2*numel(lab)));
ntrees = 100;
bagfract = bagsize / (numel(lab)-numel(valset));
[labnf, laberr] = rfSelectPoolSize(data, lab, bagfract, ntrees, ...
    maxtree, nfs, valset, w);

ntrees = 500;
[dt_lab, laboob] = rfTrainForest(data, lab, ...
    bagfract, ntrees, maxtree, labnf, [], w);

save('~/data/rfexperiments/msrcRandForestJointClassifiers.mat', ...
    'maxdata', 'maxtree', 'nfs', 'segnf', 'segerr', 'labnf', 'laberr', ...
    'dt_seg', 'segoob', 'dt_lab', 'laboob');

end


%% Testing

if DO_TEST

if 0
    load(fullfile(datadir, 'msrc_imsegs.mat'));
    load(fullfile(datadir, 'testLabels.mat'));   
    load(fullfile(datadir, 'msMultipleSegmentations_tu_tmp.mat'));
end

pg = msTestRf(imsegs(test), segfeatures(test, :), smaps(test), dt_lab, dt_seg, 0);
ignore = [5 8];  
DO_SP = 1;
[acc, cm, classcount] = msAnalyzeResult(imsegs(test), labels(test), pg, DO_SP, ignore);
cm(ignore, :) = []; cm(:, ignore) = [];        

save('~/data/rfexperiments/msrcResult_rf.mat', 'pg', 'acc', 'cm');

end