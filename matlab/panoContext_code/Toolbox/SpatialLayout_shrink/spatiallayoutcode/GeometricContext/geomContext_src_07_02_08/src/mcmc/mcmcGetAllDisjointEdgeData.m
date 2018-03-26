function [efeaturesdj, adjlistdj] = mcmcGetAllEdgeData(spfeatures, imsegs, adjlist)

for f = 1:numel(imsegs)    
    disp([num2str(f) ': getting edge data'])
    npertype = imsegs(f).nseg*4;
    [efeaturesdj{f}, adjlistdj{f}] = mcmcGetDisjointEdgeData(imsegs(f), spfeatures{f},  npertype, adjlist{f});
    disp(num2str([imsegs(f).nseg size(adjlist{f}, 1) size(adjlistdj{f},1)]))
end