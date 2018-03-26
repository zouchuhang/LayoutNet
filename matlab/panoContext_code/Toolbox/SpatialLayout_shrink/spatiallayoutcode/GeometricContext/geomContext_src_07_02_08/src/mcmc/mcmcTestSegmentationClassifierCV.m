function errors = mcmcTestSegmentationClassifierCV(classifier, data, labels, segimind, nimages)

ncv = numel(classifier);

disp([num2str(ncv) '-fold cross-validation']);

conf = zeros(size(labels));

for k = 1:ncv
    
    testind = [floor((k-1)*nimages/ncv)+1:floor(k*nimages/ncv)]; 
    keep = zeros(size(labels));
    for j = 1:numel(testind)
        keep(find(segimind==testind(j))) = 1;
    end
    ind = find(keep);       
    
    conf(ind) = test_boosted_dt_mc(classifier(k), data(ind, :));

end
    
conf = 1 ./ (1+exp(-conf));
error = mean((conf>0.5)~=(labels+1)/2);
    
disp(['error: ' num2str(error)])


ind1 = find(labels==-1);
ind2 = find(labels==1);
px = [0.025:0.05:0.975];
f1 = ksdensity(conf(ind1), px, 'support', [0 1]);
f1 = f1 / sum(f1);
f2 = ksdensity(conf(ind2), px, 'support', [0 1]);
f2 = f2 / sum(f2);
figure(1), hold off, plot(px, f1, 'r', 'LineWidth', 1);
hold on, plot(px, f2, 'g', 'LineWidth', 1);

disp(['ave conf: ' num2str(mean([1-conf(ind1) ; conf(ind2)]))])

figure(2), plot(px, f2*numel(ind2) ./ (f1*numel(ind1)+f2*numel(ind2)))
hold on, plot(px, px, '--k')
errors.err = error;
errors.pneg = f1;
errors.ppos = f2;
errors.px = px;    