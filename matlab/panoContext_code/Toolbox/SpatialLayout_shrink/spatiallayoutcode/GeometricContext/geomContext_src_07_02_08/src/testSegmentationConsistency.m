function results = testSegmentationConsistency(imsegs, smaps)

allvgoodprc = 0;
allvmajorprc = 0;
allvevergood = 0;
allvevermajor = 0;
allhgoodprc = 0;
allhmajorprc = 0;
allhevergood = 0;
allhevermajor = 0;

vtotal = 0;
htotal = 0;

nimages = numel(imsegs);

for f = 1:nimages       
       
    npixels = imsegs(f).npixels(:) / sum(imsegs(f).npixels);
    vlabels = imsegs(f).vert_labels(:);
    hlabels = imsegs(f).horz_labels(:);
    
    isevergoodv = zeros(imsegs(f).nseg, 1);
    isevermajorv = zeros(imsegs(f).nseg, 1);
    
    isevergoodh = zeros(imsegs(f).nseg, 1);
    isevermajorh = zeros(imsegs(f).nseg, 1);    
    
    ngoodv = 0;
    nmajorv = 0;

    ngoodh = 0;
    nmajorh = 0;    
    
    if size(smaps{f}, 2)==imsegs(f).nseg
        smaps{f} = smaps{f}';
    end
    
    nmaps = size(smaps{f}, 2);

  
    
    for m = 1:nmaps
        
        for s = 1:max(smaps{f}(:, m))
            vcount = zeros(1, 3);
            hcount = zeros(1, 5);
            ind = find(smaps{f}(:, m)==s);                    
            
            nv = sum((vlabels(ind)>0).*npixels(ind));
            nh = max(sum((hlabels(ind)>0).*npixels(ind)), 1E-6);              
            
            for k = 1:3                
                if ~isempty(ind)
                    vcount(k) = vcount(k) + sum((vlabels(ind)==k).*npixels(ind));
                end
            end
            for k = 1:5
                if ~isempty(ind)
                    hcount(k) = hcount(k) + sum((hlabels(ind)==k).*npixels(ind));
                end
            end                                    
            
            vcount = vcount / max(sum(vcount), 1E-10);
            hcount = hcount / max(sum(hcount), 1E-10);
        
            [maxvalv, maxlabv] = max(vcount);
            [maxvalh, maxlabh] = max(hcount);            
            
            nmajorv = nmajorv + max(vcount)*nv;
            nmajorh = nmajorh + max(hcount)*nh;
            
            ngoodv = ngoodv + (max(vcount)>0.99)*nv;
            ngoodh = ngoodh + (max(hcount)>0.99)*nh;            
            
            if max(vcount)>0.99
                isevergoodv(ind(find(vlabels(ind)>0))) = 1;
            end
            if max(hcount) > 0.99
                isevergoodh(ind(find(hlabels(ind)>0))) = 1;
            end            
            
            ind2v = find(vlabels(ind)==maxlabv);
            ind2h = find(hlabels(ind)==maxlabh);
            
            isevermajorv(ind(ind2v)) = 1;
            isevermajorh(ind(ind2h)) = 1;
        end
            
    end
    
    nv = sum((vlabels(:)>0).*npixels(:));
    nh = max(sum((hlabels(:)>0).*npixels(:)), 1E-10);  
    
    vtotal = vtotal + nv;
    htotal = htotal + nh;
    
    allvgoodprc = allvgoodprc + ngoodv/nv/nmaps/nimages;
    allhgoodprc = allhgoodprc + ngoodh/nh/nmaps/nimages;
    allvmajorprc = allvmajorprc + nmajorv/nv/nmaps/nimages;
    allhmajorprc = allhmajorprc + nmajorh/nh/nmaps/nimages;    
    allvevergood = allvevergood + sum(isevergoodv.*npixels(:));    
    allvevermajor = allvevermajor + sum(isevermajorv.*npixels(:));
    allhevergood = allhevergood + sum(isevergoodh.*npixels(:));
    allhevermajor = allhevermajor + sum(isevermajorh.*npixels(:));
    
end


allvevergood = allvevergood/vtotal;    
allvevermajor = allvevermajor/vtotal;
allhevergood = allhevergood/htotal;
allhevermajor = allhevermajor/htotal;

results.allvgoodprc = allvgoodprc;
results.allvmajorprc = allvmajorprc;
results.allvevergood = allvevergood;
results.allvevermajor = allvevermajor;
results.allhgoodprc = allhgoodprc;
results.allhmajorprc = allhmajorprc;
results.allhevergood = allhevergood;
results.allhevermajor = allhevermajor;
    