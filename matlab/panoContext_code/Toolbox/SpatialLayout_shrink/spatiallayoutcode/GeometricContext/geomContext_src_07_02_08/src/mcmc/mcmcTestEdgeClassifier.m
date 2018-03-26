function errors = mcmcTestEdgeClassifier(efeatures, adjlist, imsegs, eclassifier)

labels = {imsegs(:).labels};

[edata, elab] = formatData(efeatures, adjlist, labels);

econf = test_boosted_dt_mc(eclassifier, edata);
econf = 1 ./ (1+exp(-econf));
mean(elab)
eerror = mean((econf>0.5)~=elab);

econf2 = 1-abs(elab-econf);

ind1 = find(elab==0);
ind2 = find(elab==1);
px = [0.025:0.05:0.975];
f1 = ksdensity(econf(ind1), px, 'support', [0 1]);
f2 = ksdensity(econf(ind2), px, 'support', [0 1]);
fc = ksdensity(econf2, px, 'support', [0 1]);
%fc = fc;

errors.err = eerror;
errors.pneg = f1;
errors.ppos = f2;
errors.conf = fc;
errors.px = px;


medFS = 18;
bigFS = 20;

figure(1), hold off, plot(px, fc, 'r', 'LineWidth', 2);
%axis([0 1 0 1])
xlabel('Confidence in True Label', 'FontSize', medFS)
ylabel('Frequency', 'FontSize', medFS)
title('Same Label Confidence', 'FontSize', bigFS) 
set(gca, 'FontSize', medFS)

figure(2), hold off, plot(px, f2 ./ (f1+f2), 'LineWidth', 2)
hold on, plot(px, px, '--k')
axis([0 1 0 1])
xlabel('Estimated Probability', 'FontSize', medFS)
ylabel('Empirical Probability', 'FontSize', medFS)
%title('Same Label Confidence', 'FontSize', bigFS) 
set(gca, 'FontSize', medFS)






function [data, lab] = formatData(features, adjlist, labels)
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
    lab(c+1:c+cf) = (labels{f}(s1)==labels{f}(s2));
    c = c + cf;
end





