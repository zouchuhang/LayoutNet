function [segclassifier] = mcmcTrainSegmentClassifier2(features, labels, weights, maxdata, classparams)

if exist('classparams') && ~isempty(classparams)
    nnodes = classparams(1);
    ntrees = classparams(2);
    stopval = classparams(3);
else
    ntrees = 20;
    nnodes = 8;
    stopval = 0;
end

if ~exist('maxdata') || isempty(maxdata)
    maxdata = 25000;
end

[data, lab, w] = formatData(features, labels, weights, maxdata);

segclassifier = train_boosted_dt_2c(data, [], lab, ntrees, nnodes, stopval, w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data, lab, w] = formatData(features, labels, weights, maxdata)
% concatenate data
%disp(num2str(maxdata))
nimages = numel(features);

[tmp, nvars] = size(features{1});

% count segments
tseg = 0;
nseg = 0;
for f = 1:nimages
    nseg = nseg + sum((labels{f}~=0) & (weights{f}>0));
    tseg = tseg + numel(labels{f});
end
%disp(num2str([nseg tseg]))

data = zeros(nseg, nvars);
lab = zeros(nseg, 1);
w = zeros(nseg, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = find((labels{f}~=0) & (weights{f}>0));
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    lab(vc+1:vc+numel(ind)) = (labels{f}(ind)>0)*2-1;
    w(vc+1:vc+numel(ind)) = weights{f}(ind); % weight according to area in image
    vc = vc + numel(ind);
end

if nseg > maxdata
    rind = randperm(nseg);
    rind = rind(1:maxdata);
    data = data(rind, :);
    lab = lab(rind);
    w = w(rind);
end

w = w / sum(w);