function [vclassifier, hclassifier] = mcmcTrainSegmentClassifier2(features, labels, weights, maxdata, classparams)

if exist('classparams') && ~isempty(classparams)
    nnodes = classparams(1);
    ntrees = classparams(2);
    stopval = classparams(3);
else
    nnodes = 8;
    ntrees = 20;
    stopval = 0;
end

if ~exist('maxdata') || isempty(maxdata)
    maxdata = 25000;
end

if ~exist('weights') || isempty(weights)
    for f = 1:numel(labels)
        weights{f} = ones(size(labels{f})) / numel(labels{f});
    end
end

vlabels = labels;
hlabels = labels;
for f = 1:numel(labels)
    vlabels{f} = 1*(labels{f}==1) + 2*((labels{f}>1) & (labels{f}<7)) + ...
        3*(labels{f}==7);
    hlabels{f} = (labels{f}-1).*((labels{f}>1) & (labels{f}<7));
end

[vdata, vlab, vw] = formatData(features, vlabels, weights, maxdata);
[hdata, hlab, hw] = formatData(features, hlabels, weights, maxdata);

vnames = {'000', '090', 'sky'};
hnames = {'045', '090', '135', 'por', 'sol'};

%if size(vdata, 2) > 75
catids = []; % 75
%end
%num2str(find(std(vdata, 1)==0))
vclassifier =  train_boosted_dt_mc(vdata, catids, vnames(vlab)', ntrees, nnodes, stopval, vw, vnames);

hclassifier =  train_boosted_dt_mc(hdata, catids, hnames(hlab)', ntrees, nnodes, stopval, hw, hnames);


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