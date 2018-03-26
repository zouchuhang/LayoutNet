function [vclassifier, hclassifier] = mcmcTrainSegmentRegressor(features, labelprobs, weights)

ntrees = 20;
nnodes = 8;

maxdata = 30000;

for f = 1:numel(labelprobs)
    vprobs{f} = [labelprobs{f}(:, 1) sum(labelprobs{f}(:, 2:6), 2) labelprobs{f}(:, 7)];
    ind = find(vprobs{f}(:, 2)~=0);
    hprobs{f} = zeros(size(vprobs{f}, 1), 5); 
    hprobs{f}(ind, :) = labelprobs{f}(ind, 2:6) ./ repmat(vprobs{f}(ind, 2), [1 5]);
end

[vdata, vlab, vw] = formatData(features, vprobs, weights, maxdata);
disp(numel(vw))
[hdata, hlab, hw] = formatData(features, hprobs, weights, maxdata);
disp(numel(hw))

vnames = {'000', '090', 'sky'};
hnames = {'045', '090', '135', 'por', 'sol'};

catids = 75;

vclassifier =  train_boosted_dt_mc(vdata, catids, vnames(vlab)', ntrees, nnodes, vw, vnames);

hclassifier =  train_boosted_dt_mc(hdata, catids, hnames(hlab)', ntrees, nnodes, hw, hnames);



function [data, lab, w] = formatData(features, probs, weights, maxdata)
% concatenate data

nimages = numel(features);
nlabels = size(probs{1}, 2);

[tmp, nvars] = size(features{1});

% count segments
nseg = 0;
for f = 1:nimages
    nseg = nseg + size(features{f}, 1);
end
disp(num2str(nseg))

data = zeros(nseg, nvars);
labprob = zeros(nseg, nlabels);
w = zeros(nseg, 1);

% concatenate data
vc = 0;
for f = 1:nimages
    ind = [1:size(features{f}, 1)];
    data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
    labprob(vc+1:vc+numel(ind), :) = probs{f}(ind, :);
    w(vc+1:vc+numel(ind)) = weights{f}(ind);
    vc = vc + numel(ind);
%     for y = 1:nlabels        
%         data(vc+1:vc+numel(ind), :) = features{f}(ind, :);    
%         lab(vc+1:vc+numel(ind)) = y;
%         % weight according to (image area) * (label percentage)
%         w(vc+1:vc+numel(ind)) = weights{f}(ind).*probs{f}(ind, y); 
%         vc = vc + numel(ind);
%     end
end

if nseg > maxdata
    rind = randperm(nseg);
    data = data(rind(1:maxdata), :);
    labprob = labprob(rind(1:maxdata), :);
    w = w(rind(1:maxdata));
    nseg = maxdata;
end
tmpw = w;

w = zeros(nseg*nlabels, 1);
lab = zeros(nseg*nlabels, 1);

data = repmat(data, [nlabels 1]);
for y = 1:nlabels
    lab((y-1)*nseg+1:y*nseg) = y;
    w((y-1)*nseg+1:y*nseg) = tmpw .* labprob(:, y);
end

% remove zero weight data
ind = find(w==0);
data(ind, :) = [];
lab(ind) = [];
w(ind) = [];



w = w / sum(w);