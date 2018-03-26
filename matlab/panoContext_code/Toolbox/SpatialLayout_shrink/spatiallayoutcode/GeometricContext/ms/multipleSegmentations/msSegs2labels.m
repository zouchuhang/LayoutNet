function [mclab, mcprc, allprc, trainw] = msSegs2labels(imsegs, smap, labels, nclasses)


nseg = max(smap);

mclab = zeros(nseg, 1);
mcprc = zeros(nseg, 1);
allprc = zeros(nseg, nclasses);
trainw = zeros(nseg, 1);

imsegs.npixels = imsegs.npixels(:);

for s = 1:nseg
    ind = find(smap==s);
    for k = 1:nclasses        
        allprc(s, k) = allprc(s, k) + sum((labels(ind)==k).*imsegs.npixels(ind));        
    end    
    allprc(s, :) = allprc(s, :) / max(sum(allprc(s, :)),1E-6);
    [mcprc(s), mclab(s)] = max(allprc(s, :));
    trainw(s) = sum(imsegs.npixels(ind))/sum(imsegs.npixels);
end
