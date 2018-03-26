function [vclassifier, hclassifier] = mcmcTrainSegmentClassifier3(features, vlabels, hlabels, vw, hw)


ntrees = 20;
nnodes = 16;

maxdata = 30000;

[vdata, vlab, vw] = formatData(features, vlabels, vw, maxdata);
[hdata, hlab, hw] = formatData(features, hlabels, hw, maxdata);

vnames = {'000', '090', 'sky'};
hnames = {'045', '090', '135', 'por', 'sol'};

catids = 75;

vclassifier =  train_boosted_dt_mc(vdata, catids, vnames(vlab)', ntrees, nnodes, vw, vnames);

hclassifier =  train_boosted_dt_mc(hdata, catids, hnames(hlab)', ntrees, nnodes, hw, hnames);



function [data, lab, w] = formatData(features, labels, weights, maxdata)
% concatenate data

nimages = numel(features);

[tmp, nvars] = size(features{1});

% count segments
nseg = 0;
for f = 1:nimages
    nseg = nseg + sum(labels{f}>0);
end
disp(num2str(nseg))

data = zeros(nseg, nvars);
lab = zeros(nseg, 1);
w = zeros(nseg, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = find(labels{f}>0);
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    lab(vc+1:vc+numel(ind)) = labels{f}(ind);
    w(vc+1:vc+numel(ind)) = weights{f}(ind); % weight according to area in image
    vc = vc + numel(ind);
end
disp([vc nseg])

if nseg > maxdata
    rind = randperm(nseg);
    rind = rind(1:maxdata);
    data = data(rind, :);
    lab = lab(rind);
    w = w(rind);
    nseg = maxdata;
end

w = w / sum(w);