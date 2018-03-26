function [data, labels, imind] = mcmcGenerateRandomSegments3(imsegs, imdir, ...
            adjlist, spdata, edata, vclassifierSP, hclassifierSP, eclassifier)  
% Generate random segments for training good (1) vs. bad (-1) segments
% segment, where good segments consist entirely of one label
        
nimages = numel(imsegs);
nfeatures = 36;

nsegments = 60000;

disp(num2str(nimages))

data = zeros(nsegments, nfeatures);

labels = zeros(nsegments, 1);

imind = zeros(nsegments, 1);

count = 0;

for f = 1:nimages

    disp([num2str(f) ': ' imsegs(f).imname])
    
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    %imdata = mcmcComputeImageData(im, imsegs(f)) ;
    
    segimage = imsegs(f).segimage;
    
    %features{f} = mcmcGetSegmentFeatures(imsegs(f), spfeatures{f}, imdata,
    %gtmaps{f}, [1:max(gtmaps{f})]);
    
    [pvSP, phSP, pE, smap] = mcmcInitialize(spdata{f}, edata{f}, ...
        adjlist{f}, imsegs(f), vclassifierSP, hclassifierSP, eclassifier, 'labels');

    origmap = smap;
    
    while count < round(nsegments*(f/nimages))            
        
        %figure(1), hold off, imshow(label2rgb(origmap(segimage)))
        [smap, newseg, origseg, segadjmat] =  ...
            mcmcGenerateProposals(pE, adjlist{f}, smap, im, segimage);
        
        % display segments
        %figure(2), hold off, imagesc(label2rgb(smap(segimage)));
        %figure(1), displaySegmentGraph(im, smap(segimage), segadjmat);        
        
        neighborseg = [newseg find(segadjmat(newseg, :))];
           
        % keep all potential segments as data points
        sind = find(smap==newseg); 
        for ni = neighborseg
            smap2 = smap;
            smap2(sind) = ni;
            if ni~=newseg
                sind2 = find(smap2 > newseg);
                smap2(sind2) = smap2(sind2) - 1;
                if ni > newseg
                    ni = ni - 1;
                end
            end
            count = count + 1;
            
            data(count, :) = mcmcGetSegmentationFeatures(pvSP, phSP, pE, adjlist{f}, imsegs(f).npixels, smap2, ni);
            %data(count, :) = mcmcGetSegmentFeatures(imsegs(f), spdata{f}, imdata, smap2, ni);
            
            labels(count) = getMixUniLabel(imsegs(f), smap2, ni); % 0 for mix, 1 for uni
            
            imind(count) = f;
            
            if count == nsegments 
                break;
            end
            
        end      
        
        % remove possibility of returning to previous state
        if numel(neighborseg)>1
            neighborseg(find(neighborseg==origseg)) = [];
        end
        
        % randomly select neighbor
        nn = numel(neighborseg);
        goodind = find(labels((count-nn+1):count)); % good neighbors
        w = 0.25*ones(nn, 1);
        w(goodind) = 0.5; % so that a uniform segment is twice as likely to be picked
        w = cumsum(w / sum(w));
        ri = find(rand(1) < w);           
        
        ni = neighborseg(ri(1)); % get neighbor segment
        
        % attach segment newseg to segment ni (unless newseg==ni)
        smap(sind) = ni;
        if ni~=newseg
            sind2 = find(smap > newseg);
            smap(sind2) = smap(sind2) - 1;           
        end        
        %figure(3), hist(smap, [1:max(smap)])        

%         mask = zeros([size(im, 1) size(im, 2)]);       
%         mask(find(smap(segimage)==ni)) = 1;        
%         disp(num2str([max(smap) labels(count-nn+ri(1))]))
%         figure(2), hold off, imagesc(im .* repmat(mask, [1 1 3])), axis image                         
    end
    
    disp(num2str(mean(labels(1:count)==1)))
    
end

% remove segments that are not clearly good or bad
ind = find(labels==0);
labels(ind) = [];
data(ind, :) = [];
imind(ind) = [];
% end mcmcGenerateRandomSegments3

            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function label = getMixUniLabel(imsegs, smap, si)
% returns -1 for mixed; 0 for neither; 1 for uniform

lcount = zeros(numel(imsegs.label_names), 1);
sind = find(smap==si);

for k = sind'    
    slab = imsegs.labels(k);
    if slab > 0
        lcount(slab) = lcount(slab) + imsegs.npixels(k);
    end
end
npix = sum(lcount);

if npix > 0 
    lcount = lcount / npix;
else
    label = 0;
    return;
end

label = 0;
if max(lcount) < 0.95 || ((1-max(lcount))*npix > 500)
    label = -1;
elseif max(lcount) > 0.99
    label = 1;
end
    
    
    

        