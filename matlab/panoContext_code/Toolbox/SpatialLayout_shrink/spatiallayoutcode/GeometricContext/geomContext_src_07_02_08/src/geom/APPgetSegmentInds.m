function sinds = APPgetSegmentInds(imsegs);
% Gets pixel indices for each superpixel in imsegs

ns = imsegs.nseg;
sinds = cell(ns, 1);
for s = 1:ns
    sinds{s} = find(imsegs.segimage==s);
end
