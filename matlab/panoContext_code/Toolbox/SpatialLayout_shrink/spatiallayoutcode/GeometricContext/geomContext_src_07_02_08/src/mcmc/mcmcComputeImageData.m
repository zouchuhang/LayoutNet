function imdata = mcmcComputeImageData(im, imsegs)

[imh, imw] = size(imsegs.segimage);
imdata.yim = 1-repmat([(0:imh-1)/(imh-1)]', 1, imw);
imdata.xim = repmat([(0:imw-1)/(imw-1)], imh, 1);

[gx, gy] = gradient(im);
imdata.gradim = sqrt(sum(gx.^2 + gy.^2, 3));

minEdgeLen = sqrt(imh^2+imw^2)*0.02;
[vpdata.lines, vpdata.spinfo] = ...
    APPgetLargeConnectedEdges(rgb2gray(im), minEdgeLen, imsegs);
   [vpdata.v, vpdata.vars, vpdata.p, vpdata.hpos] = ...
       APPestimateVp(vpdata.lines, [imh imw], 0);
%  [vpdata.v vpdata.vars vpdata.p vpdata.hpos]=AppgetVP(vpdata.lines,[imh imw],0); 


imdata.vpdata = vpdata;

% get pixels in each superpixel
stats = regionprops(imsegs.segimage, 'PixelIdxList');
imdata.pixlist = {stats.PixelIdxList};     

%%by varsha%%%

for segno=1:imsegs.nseg
   BW=(imsegs.segimage==segno);
   imdata.tracedbndy{segno}= bwboundaries(BW);
end

%%%%%%%%%%%
% get boundary pixels for each superpixel
[tx, ty] = gradient(double(imsegs.segimage));
segedgeim = (tx~=0) | (ty ~=0);
stats = regionprops(double(imsegs.segimage).*segedgeim, 'PixelIdxList');
imdata.bndpixlist = {stats.PixelIdxList};
imdata.bndnpix = zeros(size(imdata.bndpixlist));
for s = 1:imsegs.nseg
    try
        imdata.bndnpix(s) = numel(imdata.bndpixlist{s});
    catch
        imdata.bndnpix(s) = 0;
    end
end

% get a random subset of pixels and boundary pixels to speed computation of
% some features
for s = 1:imsegs.nseg
    npix = imsegs.npixels(s);
    nbndpix = imdata.bndnpix(s);
    imdata.nsmpix(s) = min(npix, 1000);
    imdata.nsmbndpix(s) = floor(nbndpix/4);
    if npix>1000
        rind = randperm(npix);
        imdata.smpixlist{s} = imdata.pixlist{s}(rind(1:1000));
    else
        imdata.smpixlist(s) = imdata.pixlist(s);
    end
    rind = randperm(nbndpix);    
    try
        imdata.smbndpixlist{s} = imdata.bndpixlist{s}(rind(1:imdata.nsmbndpix(s)));
    catch
        imdata.smbndpixlist{s} = [];
    end
end