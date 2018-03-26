function [pg, pgIter] = mcmcTestImage(im, imsegs, vclassifierSP, hclassifierSP, ...
    eclassifier, vclassifier, hclassifier, segclassifier, niter, spdata, adjlist, edata)
% Computes the marginals of the geometry for the given input image
% spdata, adjlist, edata are optional inputs

DO_LABEL = 1;

grayim = rgb2gray(im);

if ~exist('spdata') || isempty(spdata)
    spdata = mcmcGetSuperpixelData(im, imsegs); 
end

if ~exist('adjlist') || ~exist('edata') || isempty(adjlist) || isempty(edata)
    [edata, adjlist] = mcmcGetEdgeData(imsegs, spdata);
end    

[pvSP, phSP, pE, smap] = mcmcInitialize(spdata, edata, ...
    adjlist, imsegs, vclassifierSP, hclassifierSP, eclassifier, 'label');


imdata = mcmcComputeImageData(im, imsegs);
imdata.pvSP = pvSP;
imdata.phSP = phSP;

nseg = max(smap);

segProb = computeSegProb(segclassifier, pvSP, phSP, pE, adjlist, smap, 1:nseg);

pgIter = zeros(imsegs.nseg, 7, niter);
pgIter(:, :, 1) = updateLabProb(vclassifier, hclassifier, imsegs, spdata, ...
    imdata, smap, 1:nseg, pgIter(:, :, 1));

figure(1), hold off, displayMarginal(mean(pgIter(:, :, 1), 3), imsegs.segimage, grayim);
%pause(1)
figure(2), hold off, displaySegments(smap, imsegs.segimage, grayim)
pause(3)

bestmap = smap;
bestprob = sum(log(segProb))+log(lognpdf(max(smap), 2.55, 0.69));

disp(num2str([max(smap) bestprob]))

% do mcmc
for iter = 2:niter                 
    
    [smap, newseg, origseg, segadjmat] = mcmcGenerateProposals(pE, adjlist, smap);
   % disp(num2str(iter))
    
    nseg2 = max(smap);
    
    tmpSegProb = segProb;    
    if nseg2 > nseg        
        tmpSegProb(newseg) = computeSegProb(segclassifier, pvSP, phSP, pE, adjlist, smap, newseg);
        if origseg~=newseg
            tmpSegProb(origseg) = computeSegProb(segclassifier, pvSP, phSP, pE, adjlist, smap, origseg);
        end
    end    
    
    % get all possible segmentation labelings for newseg
    neighborseg = [newseg find(segadjmat(newseg, :))];
    nn = numel(neighborseg);
    
    % compute transition probability to each neighbor
    transprob = zeros(nn, 1);
    
    sind = find(smap==newseg); 
    adj1 = find((smap(adjlist(:, 1))==newseg) | (smap(adjlist(:, 2))==newseg));
    for ni = 1:nn
        
        s = neighborseg(ni);        
        
        if s==newseg
            transprob(ni) = 1*lognpdf(nseg2, 2.55, 0.69);
        else
            eind = adj1(find((smap(adjlist(adj1, 1))==s) | (smap(adjlist(adj1, 2))==s)));
            smap2 = smap;
            smap2(sind) = s;
            
            sind2 = find(smap2 > newseg);
            smap2(sind2) = smap2(sind2) - 1;
            s2 = s - (s > newseg);            

            probni = computeSegProb(segclassifier, pvSP, phSP, pE, adjlist, smap2, s2);                        
            transprob(ni) = probni / tmpSegProb(newseg) / tmpSegProb(s);                        
            transprob(ni) = transprob(ni)*prod(1-pE(eind))*lognpdf(nseg2-1, 2.55, 0.69);
        end
        
    end

    %disp(num2str(transprob'))
    
    % randomly select neighbor according to transprob
    transprob = cumsum(transprob / sum(transprob));    
    s = neighborseg(find(rand(1) < transprob));
    s = s(1);
    
    sind = find(smap==newseg);
    pgIter(:, :, iter) = pgIter(:, :, iter-1);
    segProb = tmpSegProb;    
    if origseg==newseg && s==newseg  % nothing changed
        % do nothing
    elseif origseg~=newseg && s==newseg  % newseg split from origseg     
        if DO_LABEL
        pgIter(:, :, iter) = updateLabProb(vclassifier, hclassifier, ...
            imsegs, spdata, imdata, smap, [newseg origseg], pgIter(:, :, iter));  
        end
    else % s~=newseg - newseg merged into s

        smap(sind) = s;
        sind2 = find(smap > newseg);
        smap(sind2) = smap(sind2) - 1;
        s = s - (s > newseg);
        
        segProb(1:end-1) = segProb([1:newseg-1 newseg+1:end]);
        segProb(s) = computeSegProb(segclassifier, pvSP, phSP, pE, adjlist, smap, s);
        
        if DO_LABEL
        if origseg==newseg % origseg is gone, entirely merged to ni                        
            pgIter(:, :, iter) = updateLabProb(vclassifier, hclassifier, ...
                imsegs, spdata, imdata, smap, s, pgIter(:, :, iter));                                    
        else % part of origseg remains, newseg merged to ni
            pgIter(:, :, iter) = updateLabProb(vclassifier, hclassifier, ...
                imsegs, spdata, imdata, smap, [s origseg], pgIter(:, :, iter));               
        end
        end
    end
    
    currprob = sum(log(segProb))+log(lognpdf(max(smap), 2.55, 0.69));
    if mod(iter, 250)==0
        disp(num2str([iter max(smap) currprob]))        
        figure(1), displayMarginal(mean(pgIter(:, :, 1:iter), 3), imsegs.segimage, grayim);
        figure(2), hold off,  displaySegments(smap, imsegs.segimage, grayim)
        pause(0.2)
    end        
    
    if currprob > bestprob
        bestprob = currprob;
        bestmap = smap;
    end
    
end
    
    
pg = mean(pgIter, 3);    
            
disp(num2str(bestprob))
figure(4), hold off, imagesc(label2rgb(bestmap(imsegs.segimage)))
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function segProb = computeSegProb(segclassifier, pvSP, phSP, pE, adjlist, smap, sind)

nseg = numel(sind);
segProb = zeros(nseg, 1);

for i = 1:nseg
    
    s = sind(i);
    
    data = mcmcGetSegmentationFeatures(pvSP, phSP, pE, adjlist, smap, s);
    conf = test_boosted_dt_mc(segclassifier, data);
    segProb(i) = 1 / (1+exp(-conf));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function spprob = updateLabProb(vclassifier, hclassifier, imsegs, spdata, imdata, smap, sind, spprob)

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
    tmpprob = [vconf(1) vconf(2)*hconf vconf(3)];
    
    spind = find(smap==s);
    spprob(spind, :) = repmat(tmpprob, numel(spind), 1);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displaySegments(smap, segim, grayim)

pixim = smap(segim);

[gx, gy] = gradient(double(pixim));
edgepix = find(gx~=0 | gy~=0);

im = repmat(grayim, [1 1 3]);
for b = 1:3
    im((b-1)*prod(size(segim))+edgepix) = (b==1);
end
imagesc(im), axis image


function displayMarginal(pg, segim, grayim)

p000 = pg(:, 1);
p090 = sum(pg(:, 2:6), 2);
psky = pg(:, 7);

im = zeros(size(segim, 1), size(segim, 2), 3);
im(:, :, 1) = p090(segim);
im(:, :, 2) = p000(segim);
im(:, :, 3) = psky(segim);
imagesc(im.*repmat(grayim, [1 1 3])), axis image