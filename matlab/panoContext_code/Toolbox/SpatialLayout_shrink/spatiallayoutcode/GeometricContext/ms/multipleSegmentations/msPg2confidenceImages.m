function cimages = msPg2confidenceImages(imsegs, pg)

% [cimages, cnames] = msPg2confidenceImages(imsegs, pg)
% Computes confidence maps from marginal probs results
%

cimages = cell(numel(imsegs), 1);
nclasses = size(pg{1}, 2);

for f = 1:length(imsegs)
    
    cimages{f} = zeros([imsegs(f).imsize nclasses], 'single');
    
    for k = 1:size(pg{f}, 2)
        tmppg = pg{f}(:, k);
        cimages{f}(:, :, k) = tmppg(imsegs(f).segimage);
    end
             
end
    