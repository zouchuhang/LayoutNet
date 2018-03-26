function [efeatures, adjlist] = mcmcGetAllEdgeData(spfeatures, imsegs)

for f = 1:numel(imsegs)    
    disp([num2str(f) ': getting edge data'])
    [efeatures{f}, adjlist{f}] = mcmcGetEdgeData(imsegs(f), spfeatures{f});
end