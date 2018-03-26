function [coef, bias, variance] = regressMajorityPercentage(smaps, pv, ph, imsegs)
% estimate the percentage of a region occupied by the majority label, given
% some estimate of the superpixel label likelihoods

ndata = 0;
for f = 1:numel(smaps)
    for m = 1:size(smaps{f}, 2)
        ndata = ndata + max(smaps{f}(:, m));
    end
    plab{f} = [pv{f}(:, 1) repmat(pv{f}(:, 2), [1 5]).*ph{f} pv{f}(:, 3)];    
end

disp(num2str(ndata))

y = zeros(ndata, 1); % label
x = ones(ndata, 3); % statistics for regression
c = 0;
if 1 
for f = 1:numel(smaps)
    npix = imsegs(f).npixels(:);
    labels = imsegs(f).labels(:);
    
    for m = 1:size(smaps{f}, 2)
        for s = 1:max(smaps{f}(:, m))
            ind = find(smaps{f}(:, m)==s);

            mclab = 0;
            mcprc = 0;
            for k = 1:7
                prc = sum((labels(ind)==k).*npix(ind))/(sum(npix(ind).*(labels(ind)>0))+1E-10);
                if prc > mcprc
                    mcprc = prc;
                    mclab = k;
                end
            end
            if mclab > 0
                y(c+s) = mcprc;
                x(c+s, 1) = sum(plab{f}(ind, mclab).*npix(ind))/sum(npix(ind));
                [tmp, maxlab] = max(plab{f}(ind, :), [], 2);
                x(c+s, 2) = sum((maxlab==mclab).*npix(ind)) / sum(npix(ind));
                x(c+s, 3) = numel(ind);
            end
        end
        c = c + s;
    end
end
ind = find(y==0);
y(ind) = [];
x(ind, :) = [];
save '../data/tmp.mat' x y
else
    load '../data/tmp.mat'
end

thresh1 = [0.05:0.1:0.95];
tx = sum(repmat(x(:, 1), [1 numel(thresh1)])  >= repmat(thresh1, [size(x,1) 1]), 2);
hist(tx)
for k =1:numel(thresh1)
    ind = find(tx==k);
    px(k) = mean(y(ind));
    count(k) = numel(ind);
    x(ind, 1) = px(k);
end
disp(num2str(px))
disp(num2str(count))

thresh2 = [0.05:0.1:0.95];
tx = sum(repmat(x(:, 2), [1 numel(thresh2)])  >= repmat(thresh2, [size(x,1) 1]), 2);
hist(tx)
for k =1:numel(thresh2)
    ind = find(tx==k);
    px(k) = mean(y(ind));
    count(k) = numel(ind);
    x(ind, 2) = px(k);
end
disp(num2str(px))
disp(num2str(count))

thresh3 = [1 2 3 4 5 7 10 15 20 25 35 50 100 150];
tx = sum(repmat(x(:, 3), [1 numel(thresh3)])  >= repmat(thresh3, [size(x,1) 1]), 2);
hist(tx)
for k =1:numel(thresh3)
    ind = find(tx==k);
    px(k) = mean(y(ind));
    count(k) = numel(ind);
    x(ind, 3) = px(k);
end
disp(num2str(px))
disp(num2str(count))

ylog = log((y+1E-5)./(1-y+1E-5));
xlog = log((x+1E-5)./(1-x+1E-5));

coef = fminunc(@(a) objective(a, xlog, y), [1 1 1 0]'); 

%coef = robustfit(xlog, ylog, 'logistic');
%coef = regress(ylog, [ones([size(xlog, 1) 1])  xlog]);
%coef = regress(ylog, xlog);

%ylogEst = coef(1) + xlog*coef(2:end);
ylogEst = coef(end) + xlog*coef(1:end-1);
size(ylogEst)
yEst = 1./(1+exp(-ylogEst));
disp(num2str([mean(yEst) mean(y)]))
bias = mean(yEst-y);
variance = var(yEst-y);
disp(mean(abs(yEst-y)))

function val = objective(coef, xlog, y)
ylogEst = coef(end) + xlog*coef(1:end-1);
yEst = 1./(1+exp(-ylogEst));
val = sum((yEst-y).^2);
%disp([sqrt(val/numel(y)) coef'])
