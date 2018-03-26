function dt = regressMajorityPercentageDT(smaps, segdata, imsegs)
% estimate the percentage of a region occupied by the majority label, given
% some estimate of the superpixel label likelihoods

ndata = 0;
for f = 1:numel(smaps)
    for m = 1:size(smaps{f}, 2)
        ndata = ndata + max(smaps{f}(:, m));
    end
end

y = zeros(ndata, 1); % label
x = zeros(ndata, size(segdata{1}, 2)); % statistics for regression

if 1
c = 0;
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
                x(c+s, :) = segdata{f, m}(s, :);
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

n = randperm(numel(y));
y1 = y(n(1:50000));
x1 = x(n(1:50000), :);
x2 = x(n(50001:75000), :);
y2 = y(n(50001:75000));

dt = treefit(x1,y1);
dt = treeprune(dt, 'level', max(dt.prunelist)-75);
%treedisp(dt)
[c,s,n,best] = treetest(dt,'test',x2,y2,'treesize', 'se');
disp(num2str(c(best+1)))
dt = treeprune(dt,'level',best);
treedisp(dt)

