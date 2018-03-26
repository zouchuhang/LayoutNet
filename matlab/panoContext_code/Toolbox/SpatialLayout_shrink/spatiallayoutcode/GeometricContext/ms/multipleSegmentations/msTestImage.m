function [pg, data, imsegs] = msTestImage(im, imsegs, classifiers, nsegments, normalize, ...
    smaps, spdata, adjlist, edata)
%  [pg, data, imsegs] = msTestImage(im, imsegs, classifiers, nsegments, smaps,
%  spdata, adjlist, edata)
%
% Computes the marginals of the geometry for the given input image
% spdata, adjlist, edata are optional inputs
% Note: only first three arguments are required


% if imsegs is a cell, it contains the system command to make the
% segmentation image and the filename of the segmentation image

if iscell(imsegs)
    syscall = imsegs{1};
    outfn = imsegs{2};
    system(syscall);
    imsegs = processSuperpixelImage(outfn);
end

if ~exist('normalize', 'var') || isempty(normalize)
    normalize = 1;
end


labelclassifier = classifiers.labelclassifier;
segclassifier = classifiers.segclassifier;

imdata = mcmcComputeImageData(im, imsegs);

if ~exist('spdata', 'var') || isempty(spdata)
    spdata = mcmcGetSuperpixelData(im, imsegs); 
end

if ~exist('smaps', 'var') || isempty(smaps)   

    if ~exist('adjlist', 'var') || ~exist('edata', 'var') || isempty(adjlist) || isempty(edata)
        [edata, adjlist] = mcmcGetEdgeData(imsegs, spdata);
    end    
       
    eclassifier = classifiers.eclassifier; 
    ecal = classifiers.ecal;
    if isfield(classifiers, 'vclassifierSP')
        vclassifierSP = classifiers.vclassifierSP;
        hclassifierSP = classifiers.hclassifierSP;            
        [pvSP, phSP, pE] = mcmcInitialize(spdata, edata, ...
            adjlist, imsegs, vclassifierSP, hclassifierSP, eclassifier, ecal, 'none');    
    else     
        pE = test_boosted_dt_mc(eclassifier, edata);
        pE = 1 ./ (1+exp(ecal(1)*pE+ecal(2)));
    end
    smaps = generateMultipleSegmentations2(pE, adjlist, imsegs.nseg, nsegments);
end

labdata = cell(1, size(smaps, 2));
for k = 1:size(smaps, 2)
    labdata{1, k} = mcmcGetSegmentFeatures(imsegs, spdata, imdata, smaps(:, k), (1:max(smaps(:, k))));
end
    
pg = msTest(imsegs, labdata, {smaps}, labelclassifier, segclassifier, normalize);

% nsp = imsegs.nseg;  
% 
% nclasses = size(labelclassifier.wcs, 2);
% pg = zeros(nsp, nclasses);
%     
% segs = cell(nsp, 1);
%     
% for k = 1:size(smaps, 2)
% 
%     for s = 1:max(smaps(:, k))
% 
%         [segs, ind] = checksegs(segs, smaps(:, k), s);            
% 
%         if ~isempty(ind)
% 
%             labdata = mcmcGetSegmentFeatures(imsegs, spdata, imdata, smaps(:, k), s);
%             
%             vconf = test_boosted_dt_mc(vclassifier, labdata);
%             vconf = 1 ./ (1+exp(-vconf));
%             vconf = vconf / sum(vconf); 
% 
%             hconf = test_boosted_dt_mc(hclassifier, labdata);
%             hconf = 1 ./ (1+exp(-hconf));
%             hconf = hconf / sum(hconf);            
% 
%             sconf = test_boosted_dt_mc(sclassifier, labdata);
%             sconf = 1 ./ (1+exp(-sconf));           
% 
%             pgs = [vconf(1) vconf(2)*hconf vconf(3)]*sconf;
% 
%             pg(ind, :) = pg(ind, :) + repmat(pgs, numel(ind), 1);
%         end
% 
%     end
% 
% end
%         
% pg = pg ./ max(repmat(sum(pg, 2), 1, size(pg, 2)), 0.00001);    
pg = pg{1};   

data.smaps = smaps;
data.edata = edata;
data.adjlist = adjlist;
data.spdata = spdata;
data.imdata = imdata;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [segs, ind] = checksegs(segs, map, s)
% Checks whether this segment has been seen before

ind = find(map==s);

if isempty(ind)
    return;
end

% if numel(ind)==1 % already accounted for by superpixels
%     ind = [];
%     return;
% end

oldsegs = segs{ind(1)};

for k = 1:numel(oldsegs)
    if (numel(oldsegs{k})==numel(ind)) && all(oldsegs{k}==ind)
        ind = [];
        return;
    end
end

segs{ind(1)}{end+1} = ind;