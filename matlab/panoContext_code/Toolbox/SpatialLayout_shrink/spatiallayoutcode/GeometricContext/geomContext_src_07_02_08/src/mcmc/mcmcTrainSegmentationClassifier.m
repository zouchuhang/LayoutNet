function segclassifier = mcmcTrainSegmentationClassifier(data, labels, segimind, trainind)

ntrees = 10;
nnodes = 8;

disp(num2str(mean(labels==1)))

keep = zeros(numel(labels), 1);

for k = 1:numel(trainind)
    keep(find(segimind==trainind(k))) = 1;
end
ind = find(keep);

segclassifier = train_boosted_dt_2c(data(ind, :), [], labels(ind), ntrees, nnodes);

    

