function writeConfidenceImages(imsegs, pg, outdir)

[cimages, cnames] = pg2confidenceImages(imsegs, pg);

for f = 1:numel(imsegs)
    bn = strtok(imsegs(f).imname, '.');
    for k = 1:size(cimages{f}, 3)        
        imwrite(cimages{f}(:, :, k), [outdir '/' bn '_c_' cnames{k} '.jpg'], 'Quality', 100);        
    end
end