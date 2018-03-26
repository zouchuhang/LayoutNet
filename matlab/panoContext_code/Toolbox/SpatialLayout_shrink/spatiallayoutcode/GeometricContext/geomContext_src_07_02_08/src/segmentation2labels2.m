function [vlab, hlab, vw, hw] = segmentation2labels2(imsegs, smaps)
% [vlab, hlab, trainw] = segmentation2labels2(imsegs, smaps)
% imsegs(nimages)
% smaps{nimages}(nsp, nmaps)
% vlab{nimages, nmaps}(nsp)
% hlab{nimages, nmaps}(nsp)
% vw{nimages, nmaps}(nsp)
% hw{nimages, nmaps}(nsp)

if ~iscell(smaps)
    smaps = {smaps};
end

nimages = numel(imsegs);
nmaps = size(smaps{1}, 2);

vlab = cell(nimages, nmaps);
hlab = cell(nimages, nmaps);
trainw = cell(nimages, nmaps);

for f = 1:numel(imsegs)   
    
    for m = 1:size(smaps{f}, 2)

        smap = smaps{f}(:, m);    
        nseg = max(smap);

        npixels = imsegs(f).npixels(:) / sum(imsegs(f).npixels(:));
        vlabels = imsegs(f).vert_labels(:);
        hlabels = imsegs(f).horz_labels(:);

        vlab{f,m} = zeros(nseg, 1);
        hlab{f,m} = zeros(nseg, 1);
        vw{f,m} = zeros(nseg, 1);
        hw{f,m} = zeros(nseg, 1);         
        
        for s = 1:nseg                      
            
            vcount = zeros(1, 3);
            hcount = zeros(1, 5);    
            ind = find(smap==s);
            for k = 1:3        
                vcount(k) = sum((vlabels(ind)==k).*npixels(ind));                       
            end    
            if (sum(vcount)>0) && (max(vcount/sum(vcount)) > 0.98)
                [tmp, vlab{f, m}(s)] = max(vcount);
            elseif (sum(vcount)>0) && (max(vcount/sum(vcount)) < 0.95)
                vlab{f, m}(s) = -1;
            end
            for k = 1:5        
                hcount(k) = sum((hlabels(ind)==k).*npixels(ind));                       
            end    
            
            if (sum(hcount)>0) && (max(hcount/sum(hcount)) > 0.99)
                [tmp, hlab{f, m}(s)] = max(hcount);
            elseif (sum(hcount)>0) && (max(hcount/sum(hcount)) < 0.95)
                hlab{f, m}(s) = -1;
            end           
            vw{f, m}(s) = sum(vcount);
            hw{f, m}(s) = sum(hcount);            
        end
    end
end