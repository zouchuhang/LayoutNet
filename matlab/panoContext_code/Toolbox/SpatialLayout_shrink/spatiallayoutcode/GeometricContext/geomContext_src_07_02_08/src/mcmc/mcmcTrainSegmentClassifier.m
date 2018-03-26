function [vclassifier, hclassifier] = ...
    mcmcTrainSegmentClassifier(imsegs, features, labels)

ntrees = 15;
nnodes = 8;

maxdata = 25000;

vlabels = labels;
hlabels = labels;
for f = 1:numel(labels)
    vlabels{f} = 1*(labels{f}==1) + 2*((labels{f}>1) & (labels{f}<7)) + ...
        3*(labels{f}==7);
    hlabels{f} = (labels{f}-1).*((labels{f}>1) & (labels{f}<7));
end

[vdata, vlab, vw] = formatData(features, vlabels, maxdata);
[hdata, hlab, hw] = formatData(features, hlabels, maxdata);

vnames = {'000', '090', 'sky'};
hnames = {'045', '090', '135', 'por', 'sol'};

catids = 75;

vclassifier =  train_boosted_dt_mc(vdata, catids, vnames(vlab)', ntrees, nnodes, vw, vnames);

hclassifier =  train_boosted_dt_mc(hdata, catids, hnames(hlab)', ntrees, nnodes, hw, hnames);



function [data, lab, w] = formatData(features, labels, maxdata)
% concatenate data

disp('assuming area at 54')

nimages = numel(features);

[tmp, nvars] = size(features{1});

% count segments
nseg = 0;
for f = 1:nimages
    nseg = nseg + sum(labels{f}~=0);
end

data = zeros(nseg, nvars);
lab = zeros(nseg, 1);
w = zeros(nseg, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = find(labels{f}~=0);
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    lab(vc+1:vc+numel(ind)) = labels{f}(ind);
    w(vc+1:vc+numel(ind)) = features{f}(ind, 54); % weight according to area in image
    vc = vc + numel(ind);
end

% remove segments that are very small (may be labeling error)
[f, x] = ksdensity(w, 'support', [0 1]);
figure(1), hold off, plot(x, f);
% ind = find(w < 0.001);
% disp(num2str([numel(ind) numel(w)]))
if numel(w) > maxdata
    ind = randperm(numel(w));
    ind = ind(1:maxdata);
    w = w(ind);
    data = data(ind, :);
    lab = lab(ind);
end

w = w / sum(w);