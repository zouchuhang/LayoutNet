function [vclassifier, hclassifier] = ...
    mcmcTrainSuperpixelClassifier(spfeatures, imsegs)

ntrees = 25;
nnodes = 8;
ndata = 10000;
vnames = imsegs(1).vert_names;
hnames = imsegs(1).horz_names;
spvlabels = {imsegs(:).vert_labels};
sphlabels = {imsegs(:).horz_labels};

% train {ground, vertical, sky} classifier
[vdata, vlab] = formatData(spfeatures, spvlabels, ndata); 
vclassifier = train_boosted_dt_mc(vdata, [], vnames(vlab)', ...
    ntrees, nnodes, [], vnames);

% train vertical subclass classifier  
[hdata, hlab] = formatData(spfeatures, sphlabels, ndata); 
hclassifier = train_boosted_dt_mc(hdata, [], hnames(hlab)', ...
    ntrees, nnodes, [], hnames);
    

function [data, lab] = formatData(features, labels, ndata)
% concatenate data and select ndata random datapoints

nimages = numel(features);

[tmp, nvars] = size(features{1});

% count superpixels
nsp = 0;
for f = 1:nimages
    nsp = nsp + sum(labels{f}~=0);
end

data = zeros(nsp, nvars);
lab = zeros(nsp, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = find(labels{f}~=0);
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    lab(vc+1:vc+numel(ind)) = labels{f}(ind);
    vc = vc + numel(ind);
end

% select random ndata points
if nsp > ndata
    rind = randperm(nsp);
    rind = rind(1:ndata);
    data = data(rind, :);
    lab = lab(rind);
end
