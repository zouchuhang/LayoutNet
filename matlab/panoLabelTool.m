% shell script for labeling panos
clear;
clc;
close all;

% path
addpath(genpath('./minFunc_2012/'));
addpath(genpath('./panoContext_code/'));

% data
% train & test split
config.folderName = '../panoContext_data/bedroom/';
config.annotationfile = [config.folderName 'ANNO_ALL.mat'];
load(config.annotationfile);
load([config.folderName 'IMGLIST.mat']); % image name
config.IMGLIST = IMGLIST;
clear IMGLIST
anno_valid = false( 1, length(ANNO_ALL));
for i = 1:length(anno_valid)
    anno_valid(i) = ~isempty(ANNO_ALL(i).ANNO3D) && ANNO_ALL(i).ANNO3D.b_singleroom;
end
%clear ANNO_ALL;
config.anno_valid = anno_valid;

config.valid_anno_id = find(anno_valid);
config.valid_data_id = 1:sum(anno_valid);

clear anno_valid
clear i

config.TEST_DATA_IDS = [4 5 8 14 16 18 22 25 28 33 36 38 44 48 68 72 99 111 112 118 ...
    122 124 137 142 159 196 199 201 205 211 215 222 224 231 238 255 ...
    256 264 274 275 286 315 328 333 341 342 354 363 368 370 401 404 417];
config.TRAIN_DATA_IDS = setdiff(config.valid_anno_id, config.TEST_DATA_IDS);


% options
options = [];
options.display = 'none';
options.MaxIter = 100;
options.Method = 'lbfgs';
options.LS_init = 2;

% run
for i = 1:numel(config.TEST_DATA_IDS)
    
    id = config.TEST_DATA_IDS(i);
    PanoLabelIm(id, ANNO_ALL, config, options);
end
