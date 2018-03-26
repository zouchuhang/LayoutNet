function [cimages, cnames] = pg2confidenceImages(imsegs, pg)
% [cimages, cnames] = pg2confidenceImages(imsegs, pg)
% Computes confidence maps from marginal probs results
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

cimages = cell(length(imsegs), 1);

cnames = {'000', '090-045', '090-090', '090-135', '090-por', '090-sol', 'sky', '090'};

for f = 1:length(imsegs)

    [pv, ph] = splitpg(pg{f});
    
    cimages{f} = zeros([imsegs(f).imsize 8]);
    
    for k = 1:size(pg{f}, 2)
        tmppg = pg{f}(:, k);
        cimages{f}(:, :, k) = tmppg(imsegs(f).segimage);
    end
    tmppg = pv(:, 2);
    cimages{f}(:, :, end) = tmppg(imsegs(f).segimage);
             
end
    