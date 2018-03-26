function classifier = msTrainLabelClassifier(features, labels, weights, classnames, maxdata, classparams)

if exist('classparams', 'var') && ~isempty(classparams)
    nnodes = classparams(1);
    ntrees = classparams(2);
    stopval = classparams(3);
else
    nnodes = 8;
    ntrees = 20;
    stopval = 0;
end

if ~exist('maxdata', 'var') || isempty(maxdata)
    maxdata = 25000;
end

if ~exist('weights', 'var') || isempty(weights)
    for f = 1:numel(labels)
        weights{f} = ones(size(labels{f})) / numel(labels{f});
    end
end

[data, lab, w] = formatData(features, labels, weights, maxdata);

if ~exist('classnames', 'var')
    classnames = [];
end

%if size(vdata, 2) > 75
catids = []; % 75
%end
%num2str(find(std(vdata, 1)==0))
classifier =  train_boosted_dt_mc(data, catids, lab, ntrees, nnodes, stopval, w);


%% Reformat the input data to be used by classifier 
function [data, lab, w] = formatData(features, labels, weights, maxdata)
% concatenate data

nimages = numel(features);

[tmp, nvars] = size(features{1});

% count segments
nseg = 0;
for f = 1:nimages
    nseg = nseg + sum(labels{f}>0);
end
%disp(num2str(nseg))

data = zeros(nseg, nvars);
lab = zeros(nseg, 1);
w = zeros(nseg, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = find(labels{f}>0);
    %disp(num2str([f size(labels{f}) size(features{f})]))
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    lab(vc+1:vc+numel(ind)) = labels{f}(ind);
    w(vc+1:vc+numel(ind)) = weights{f}(ind); % weight according to area in image
    vc = vc + numel(ind);
end
%disp([vc nseg])

if nseg > maxdata
    rind = randperm(nseg);
    rind = rind(1:maxdata);
    data = data(rind, :);
    lab = lab(rind);
    w = w(rind);
    nseg = maxdata;
end

w = w / sum(w);