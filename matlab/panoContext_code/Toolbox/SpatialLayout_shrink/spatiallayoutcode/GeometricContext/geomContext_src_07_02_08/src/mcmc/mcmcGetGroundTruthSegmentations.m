function [gtmaps, gtlabels] = mcmcGetGroundTruthSegmentations(imsegs, adjlist)

nregions = [];
for f = 1:numel(imsegs)
    %disp(num2str(f))
    nseg = imsegs(f).nseg;
    adjmat = zeros(size(imsegs(f).adjmat));    
    
    for k = 1:size(adjlist{f},1)   
        s1 = adjlist{f}(k, 1);
        s2 = adjlist{f}(k, 2);
        if imsegs(f).labels(s1)==imsegs(f).labels(s2)
            adjmat(s1, s2) = 1;
            adjmat(s2, s1) = 1;
        end
    end
    
    gtmaps{f} = graphComponents(adjmat);
    
    nregions(f) = max(gtmaps{f});
    
    for k = 1:max(gtmaps{f})
        ind = find(gtmaps{f}==k);
        gtlabels{f}(k) = imsegs(f).labels(ind(1));
    end
end

%hist(nregions, [1:max(nregions)]);