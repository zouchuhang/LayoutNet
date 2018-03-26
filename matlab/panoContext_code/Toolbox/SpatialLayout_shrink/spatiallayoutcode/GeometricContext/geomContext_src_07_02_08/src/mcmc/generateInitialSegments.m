function [labdata, segdata, labels, weights, smap] = generateInitialSegments(imsegs, imdir, ...
            adjlist, spdata, vclassifierSP, hclassifierSP, edata, eclassifier)  
      
ngoodv = 0;
ngoodh = 0;
ntotalv = 0;
ntotalh = 0;

for f = 1:numel(imsegs)
    disp([num2str(f) ': ' imsegs(f).imname])

    im = im2double(imread([imdir '/' imsegs(f).imname]));    
    [labdata{f}, segdata{f}, labels{f}, weights{f}, smap{f}, maxvalv, maxvalh]  = ...
        generateInitialSegmentsImage(im, imsegs(f), adjlist{f}, spdata{f}, ...
        vclassifierSP, hclassifierSP, edata{f}, eclassifier);
    ngoodv = ngoodv + sum(maxvalv.*weights{f});
    ngoodh = ngoodh + sum(maxvalh.*weights{f});
    ntotalv = ntotalv + sum(imsegs(f).npixels(:).*(imsegs(f).vert_labels(:)~=0))/sum(imsegs(f).npixels);
    ntotalh = ntotalh + sum(imsegs(f).npixels(:).*(imsegs(f).horz_labels(:)~=0))/sum(imsegs(f).npixels);
    disp(num2str([ngoodv/ntotalv ngoodh/ntotalh]));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [labdata, segdata, labels, weights, smap, maxvalv, maxvalh] = ...
    generateInitialSegmentsImage(im, imsegs, adjlist, spdata, vclassifierSP, hclassifierSP, edata, eclassifier);

[pvSP, phSP, pE, smap] = mcmcInitialize(spdata, edata, ...
    adjlist, imsegs, vclassifierSP, hclassifierSP, eclassifier, 'labels');

imdata = mcmcComputeImageData(im, imsegs);
imdata.pvSP = pvSP;
imdata.phSP = phSP;

nseg = max(smap);

labdata = mcmcGetSegmentFeatures(imsegs, spdata, imdata, smap, [1:nseg]);
segdata = mcmcGetSegmentationFeatures(pvSP, phSP, pE, adjlist, imsegs.npixels, smap, [1:nseg]);
[labels, weights, maxvalv, maxvalh] = getLabels(imsegs, smap);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [labels, weights, maxvalv, maxvalh] = getLabels(imsegs, smap)
% labels: {-1 0 1 2 3 4 5 6 7} --> {mixed, borderline, ground, left,
% center, right, porous, solid}

labels = zeros(max(smap), 1);
weights = zeros(max(smap), 1);
maxvalv = zeros(max(smap), 1);
maxvalh = zeros(max(smap), 1);

for s = 1:max(smap)

    lcount = zeros(numel(imsegs.label_names), 1);
    sind = find(smap==s);
    for k = sind   
        slab = imsegs.labels(k);
        if slab > 0
            lcount(slab) = lcount(slab) + imsegs.npixels(k);
        end
    end
    npix = sum(lcount);

    labels(s) = 0;
    weights(s) = 0;
    if npix > 0 
        lcount = lcount / npix;
        vcount = [lcount(1) sum(lcount(2:6)) lcount(7)];        
        hcount = lcount(2:6);
        [maxvalv(s), maxlabv] = max(vcount);
        [maxvalh(s), maxlabh] = max(hcount);
        if maxlabv==2
            labels(s) = 1 + maxlabh;
        else
            labels(s) = 1*(maxlabv==1) + 7*(maxlabv==3);
        end
                            
%        [maxval, maxlab] = max(lcount);                   
%        labels(s) = maxlab;
%         if maxval < 0.95 || ((1-maxval)*npix > 500)
%             labels(s) = -1;
%         elseif maxval >= 0.99
%             labels(s) = maxlab;
%         end
        weights(s) = npix / sum(imsegs.npixels);
        %maxvalv(s) = max([lcount(1) sum(lcount(2:6)) lcount(7)]);        
        %maxvalv(s) = maxvalv(s)*(maxvalv(s)>=0.95);
        if sum(lcount(2:6))>0    
            maxvalh(s) = max(lcount(2:6));        
            %maxvalh(s) = maxvalh(s)*(max(lcount(2:6))>=0.95*sum(lcount(2:6)));
        end
            
    end
    
end


        