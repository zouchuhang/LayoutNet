function [pv, ph, smap, vacc, hacc] = testImageBottomUp(im, imsegs, vclassifierSP, hclassifierSP, ...
    eclassifier, vclassifier, hclassifier, segclassifier, spdata, adjlist, edata)

% Initialize
nsp = imsegs.nseg;

grayim = rgb2gray(im);

if ~exist('spdata') || isempty(spdata)
    spdata = mcmcGetSuperpixelData(im, imsegs); 
end

if ~exist('adjlist') || ~exist('edata') || isempty(adjlist) || isempty(edata)
    [edata, adjlist] = mcmcGetEdgeData(imsegs, spdata);
end    

t = 1;

[pvSP, phSP, pE, smap{t}] = mcmcInitialize(spdata, edata, ...
    adjlist, imsegs, vclassifierSP, hclassifierSP, eclassifier, 'label');

imdata = mcmcComputeImageData(im, imsegs);
imdata.pvSP = pvSP;
imdata.phSP = phSP;

adjmat = zeros(imsegs.nseg, imsegs.nseg);

pv(:, :, t) = pvSP;
ph(:, :, t) = phSP;

while (1)

    [vmaxval(:, t), vmax(:, t)] = max(pv(:, :, t), [], 2);
    [hmaxval(:, t), hmax(:, t)] = max(ph(:, :, t), [], 2);

    ind = find(vmax{t}==2);
    if t==1 || (any(vmax(:, t)~=vmax(:, t-1)) || any(hmax(:, t)(ind)~=hmax(:, t-1)(ind)))
            
        for k = 1:size(adjlist,1)
            s1 = adjlist(k, 1);
            s2 = adjlist(k, 2);
            if (vmax(:, t)(s1)==vmax(:, t)(s2)) && (vmax(:, t)(s1)~=2 || (hmax(:, t)(s1)==hmax(:, t)(s2)))
                adjmat(s1, s2) = 1;
                adjmat(s2, s1) = 1;
            else
                adjmat(s1, s2) = 0;
                adjmat(s2, s1) = 0;
            end
        end
    
        t = t + 1;
        
        smap{t} = graphComponents(adjmat);
        
        pv(:, :, t) = pv(:, :, t-1);
        ph(:, :, t) = ph(:, :, t-1);
        if t==2
            [pv(:, :, t), ph(:, :, t)] = updateLabProb(vclassifier, hclassifier, imsegs, ...
                spdata, imdata, smap(:, t), [1:max((:, t))], pv(:, :, t), ph(:, :, t));
        else
            for seg1 = 1:max(smap(:, t))
                ind1 = find(smap(:, t)==seg1);
                seg2 = smap(:, t-1)(ind1(1));
                ind2 = find(smap(:, t-1)==seg2);
                if numel(union(ind1, ind2))~=numel(ind1)
                    [pv(:, :, t)  ph(:, :, t)] = updateLabProb(vclassifier, hclassifier, ...
                        imsegs, spdata, imdata, smap(:, t), seg1, pv(:, :, t), ph(:, :, t));
                end
            end
        end
    else
        break;
    end
    
end % end while

vlab = imsegs.vert_labels(:);
hlab = imsegs.horz_labels(:);
npixels = imsegs.npixels;

for t = 1:numel(pv)

    vacc(t) = sum((vmax(:, t)==vlab).*npixels)/sum((vlab~=0).*npixels);
    hacc(t) = sum((hmax(:, t)==hlab).*npixels)/sum((hlab~=0).*npixels);

    if t ~= numel(pv)-1
        figure((t-1)*2+1), displaySegments(smap(:, t), imsegs.segimage, grayim);
        figure((t-1)*2+2), displayMarginal(pv(:, t), imsegs.segimage, ...
            grayim*0.5+0.5*ones(size(imsegs.segimage)));
    end

end

disp(['vacc: ' num2str(vacc)])
disp(['hacc: ' num2str(hacc)])   
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pv,  ph] = updateLabProb(vclassifier, hclassifier, imsegs, spdata, imdata, smap, sind, pv, ph)

nseg = numel(sind);

for i = 1:nseg

    s = sind(i);    
    
    data = mcmcGetSegmentFeatures(imsegs, spdata, imdata, smap, s);
    vconf = test_boosted_dt_mc(vclassifier, data);
    vconf = 1 ./ (1+exp(-vconf));
    vconf = vconf / sum(vconf);    
    hconf = test_boosted_dt_mc(hclassifier, data);
    hconf = 1 ./ (1+exp(-hconf));
    hconf = hconf / sum(hconf);     
        
    spind = find(smap==s);
    pv(spind, :) = repmat(vconf, numel(spind), 1);
    ph(spind, :) = repmat(hconf, numel(spind), 1);
    
end        
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function displaySegments(smap, segim, grayim)

pixim = smap(segim);

[gx, gy] = gradient(double(pixim));
edgepix = find(gx~=0 | gy~=0);

im = repmat(grayim, [1 1 3]);
for b = 1:3
    im((b-1)*prod(size(segim))+edgepix) = (b==1);
end
imagesc(im), axis image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayMarginal(pv, segim, grayim)

p000 = pv(:, 1);
p090 = pv(:, 2);
psky = pv(:, 3);

im = zeros(size(segim, 1), size(segim, 2), 3);
im(:, :, 1) = p090(segim);
im(:, :, 2) = p000(segim);
im(:, :, 3) = psky(segim);
%imagesc(im.*repmat(grayim, [1 1 3])), axis image    
imagesc(im), axis image 