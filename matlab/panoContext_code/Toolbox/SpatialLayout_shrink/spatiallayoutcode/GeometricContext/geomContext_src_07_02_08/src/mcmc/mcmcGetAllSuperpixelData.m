function spfeatures = mcmcGetAllSuperpixelData(imdir, imsegs)

for f = 1:numel(imsegs)    
    disp([num2str(f) ': getting superpixel features'])
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    spfeatures{f} = mcmcGetSuperpixelData(im, imsegs(f)); 
end