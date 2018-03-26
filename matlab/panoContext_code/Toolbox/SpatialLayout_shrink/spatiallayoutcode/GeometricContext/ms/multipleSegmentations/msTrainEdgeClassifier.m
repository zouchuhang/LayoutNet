function eclassifier = mcmcTrainEdgeClassifier(efeatures, adjlist, imsegs, labels)

ntrees = 20;
nnodes = 8;
ndata = 50000;

% train {ground, vertical, sky} classifier
[edata, elab] = formatData(efeatures, adjlist, labels, ndata); 
mean(elab==1)
eclassifier = train_boosted_dt_2c(edata, [], elab, ntrees, nnodes, 0);




function [data, lab] = formatData(features, adjlist, labels, ndata)
% concatenate data and select ndata random datapoints

nimages = numel(features);

[tmp, nvars] = size(features{1});

% count edges
ne = 0;
for f = 1:nimages    
    ne = ne + size(features{f}, 1);
end

data = zeros(ne, nvars);
lab = zeros(ne, 1);

% concatenate data
c = 0;
for f = 1:nimages
    cf = size(features{f}, 1);    
    data(c+1:c+cf, :) = features{f};    
    s1 = adjlist{f}(:, 1);
    s2 = adjlist{f}(:, 2);
    lab(c+1:c+cf) = (labels{f}(s1)==labels{f}(s2))*2-1;
    c = c + cf;
end

% select random ndata points
if ne > ndata
    rind = randperm(ne);
    rind = rind(1:ndata);
    data = data(rind, :);
    lab = lab(rind);
end
