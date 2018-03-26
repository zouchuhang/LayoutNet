function [eparams, errors] = msCalibrateEdgeClassifier(efeatures, adjlist, imsegs, eclassifier, labels, ncv)
% [eparams, errors] = calibrateEdgeClassifier(efeatures, adjlist, imsegs,
% eclassifier, ncv)

nimages = numel(imsegs);
for k = 1:ncv    
    if ncv > 1
        testind = [(k-1)*nimages/ncv+1:k*nimages/ncv];
        trainind = setdiff([1:nimages], testind);
    else
        trainind = (1:nimages);
    end
    [edata{k}, elab{k}] = formatData(efeatures(trainind), adjlist(trainind), labels(trainind));
    econf{k} = test_boosted_dt_mc(eclassifier, edata{k});
%    econf{k} = 1 ./ (1+exp(-econf{k}));
end
   
for k = 1:ncv
    disp(['iter: ' num2str(k)])
    if ncv>1
        traink = setdiff([1:ncv], k);
    else
        traink = k;
    end
    eparams{k} = fminunc(@(x) objective(x, cat(1, econf{traink}), cat(1, elab{traink})), [-1 0], optimset('TolFun', 0.001)); 
end

for k = 1:ncv
    econf{k} = 1 ./ (1+exp(eparams{k}(1)*econf{k}+eparams{k}(2)));
end

elab = cat(1, elab{:});
econf = cat(1, econf{:});

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

figure(1), hold on, plot(px, fc, 'y', 'LineWidth', 2);
%axis([0 1 0 1])
xlabel('Confidence in True Label', 'FontSize', medFS)
ylabel('Frequency', 'FontSize', medFS)
title('Same Label Confidence', 'FontSize', bigFS) 
set(gca, 'FontSize', medFS)

figure(2), hold on, plot(px, f2 ./ (f1+f2), 'y', 'LineWidth', 2)
hold on, plot(px, px, '--k')
axis([0 1 0 1])
xlabel('Estimated Probability', 'FontSize', medFS)
ylabel('Empirical Probability', 'FontSize', medFS)
%title('Same Label Confidence', 'FontSize', bigFS) 
set(gca, 'FontSize', medFS)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = objective(param, econf, elab)
econf = 1./(1+exp(param(1)*econf+param(2)));
px = [0.025:0.05:0.975];
f1 = ksdensity(econf(elab==0), px, 'support', [0 1])+1E-10;
f2 = ksdensity(econf(elab==1), px, 'support', [0 1])+1E-10;
f1 = f1 / sum(f1+f2);
f2 = f2 / sum(f1+f2);
err = sum((f1+f2).*(px - f2./(f1+f2)).^2);
disp(num2str([sum((f1+f2).*abs(px - f2./(f1+f2))) param]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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