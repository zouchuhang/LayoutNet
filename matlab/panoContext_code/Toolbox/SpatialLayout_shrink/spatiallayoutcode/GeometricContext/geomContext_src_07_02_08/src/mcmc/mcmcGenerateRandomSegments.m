function [uniseg, mixseg] = mcmcGenerateRandomSegments(imsegs, ...
    spdata, edata, adjlist, vclassifierSP, hclassifierSP, eclassifier)

nimages = numel(imsegs);

nmix = 10000;
nuni = 10000;

cuni = 0;
cmix = 0;

% probability that any segment is selected (gives expected # of iterations)
pchoose = ((nmix+nuni)/nimages) / 1000; 

uniseg = cell(nuni, 1);
mixseg = cell(nmix, 1);

for f = 1:nimages

    disp(f)
    im = imresize(im2double(imread(['../images/all_images/' imsegs(f).imname])), 0.5);
    segimage = imresize(imsegs(f).segimage, 0.5);
    
    
    targetUni = ceil(f/nimages*nuni);
    targetMix = ceil(f/nimages*nmix);
    
    truelab = imsegs(f).labels;
    
    [pvSP, phSP, pE, smap, lmap, adjmat] = mcmcInitialize(spdata{f}, edata{f}, ...
        adjlist{f}, imsegs(f), vclassifierSP, hclassifierSP, eclassifier);
    %pvhSP = [pvSP(:, 1) repmat(pvSP(:, 2), 1, 5).*phSP pvSP(:, 3)];
    origmap = smap;
    
%     figure(1)
%     for k = 1:max(smap)
%         imshow(smap(segimage)==k)
%         pause(0.25)
%     end
    
    while (cuni < targetUni) || (cmix < targetMix)
                        
        %figure(1), hold off, imshow(label2rgb(origmap(segimage)))
        [smap, newseg, origseg, adjmat, segadjmat] =  ...
            mcmcGenerateProposals(pE, adjlist{f}, smap, adjmat);
        
        neighborseg = [newseg find(segadjmat(newseg, :))];
        rs = neighborseg(ceil(rand(1)*numel(neighborseg)));        
        
        % merge newseg with rs        
        if rs ~= newseg
            smap(find(smap==newseg)) = rs;            
            ind = find(smap>newseg);
            smap(ind) = smap(ind) - 1;            
            for k = 1:size(adjlist{f}, 1)
                s1 = adjlist{f}(k, 1);
                s2 = adjlist{f}(k, 2);
                if smap(s1)==smap(s2) 
                    adjmat(s1, s2) = 1;
                    adjmat(s2, s1) = 1;
                end
            end
            
        end       
        
        if cuni < targetUni && rand(1) < pchoose
            rs = smap(ceil(rand(1)*max(smap)));
            
            ind = find(smap==rs);
            if all(truelab(ind)==truelab(ind(1)))
                cuni = cuni + 1;
                uniseg{cuni} = ind;     
                
                mask = zeros([size(im, 1) size(im, 2)]);
                for k = 1:numel(ind)
                    mask(find(segimage==ind(k))) = 1;
                end
                figure(1), hold off, imshow(im .* repmat(mask, [1 1 3]))
                
            end
        end
       
        if cmix < targetMix && rand(1) < pchoose
            rs = smap(ceil(rand(1)*max(smap)));
            
            ind = find(smap==rs);
            if ~all(truelab(ind)==truelab(ind(1)))
                cmix = cmix + 1;
                uniseg{cmix} = ind;
                
                mask = zeros([size(im, 1) size(im, 2)]);
                for k = 1:numel(ind)
                    mask(find(segimage==ind(k))) = 1;
                end              
                figure(2), hold off, imshow(im .* repmat(mask, [1 1 3]))                     
            end                                  
            
        end        
        
    end
    
end
            
            
    
    
    

