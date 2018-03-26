function pg = msTest(imsegs, data, maps, labelclassifier, segclassifier, normalize)
% [vacc, hacc, vcm, hcm] = testMultipleSegmentationsCV2(imsegs, labdata, segdata, maps, vclassifier, hclassifier, sclassifier, ncv)
   
if ~exist('normalize', 'var') || isempty(normalize)
    normalize = 1;
end

pg = cell(numel(imsegs), 1);
allcount = 0;

nclasses = size(labelclassifier.wcs,2);

for f = 1:numel(imsegs)
            
    %disp(num2str(f))
    
    nsp = imsegs(f).nseg;   

    pg{f} = zeros(nsp, nclasses);
    
    
    segs = cell(nsp, 1);

    count = zeros(nsp, 1);
    
    if ~normalize
        sum_s = zeros(nsp, 1);
    end
    
    

    for k = 1:size(data, 2)
        if numel(data{f,k})>0
        
        isvalid = false(size(data{f, k}, 1), 1);
        for s = 1:size(data{f, k}, 1)
            [segs, ind] = checksegs(segs, maps{f}(:, k), s);            
            isvalid(s) = ~isempty(ind);
        end
            
        yconf_all = test_boosted_dt_mc(labelclassifier, data{f, k}(isvalid, :));
        if exist('segclassifier','var')
        sconf_all = test_boosted_dt_mc(segclassifier, data{f, k}(isvalid, :)); 
        end          
        
        isvalid = find(isvalid);
        for s = 1:numel(isvalid)     
            
            ind = find(maps{f}(:, k)==isvalid(s));

            
            count(ind) = count(ind) + 1;

            yconf = yconf_all(s, :);
            yconf = 1 ./ (1+exp(-yconf));
            
            if normalize
                yconf = yconf / sum(yconf);        
            end

            if exist('segclassifier','var')
            sconf = sconf_all(s, :);
            sconf = 1 ./ (1+exp(-sconf));           
            pgs = yconf*sconf;
            else
              pgs = yconf;  
            end
            
            if ~normalize
                sum_s(ind) = sum_s(ind) + sconf;
            end
            pg{f}(ind, :) = pg{f}(ind, :) + repmat(pgs, numel(ind), 1);
            
        end
          
        end
    end
     
    if normalize
        pg{f} = pg{f} ./ max(repmat(sum(pg{f}, 2), 1, size(pg{f}, 2)), 0.00001);    
    else
        pg{f} = pg{f} ./ repmat(sum_s, [1 size(pg{f}, 2)]);
    end 
        
    allcount = allcount + mean(count);
    
end
        
allcount = allcount / numel(imsegs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [segs, ind] = checksegs(segs, map, s)
% Checks whether this segment has been seen before

ind = find(map==s);

if 0

if numel(ind)==1 % already accounted for by superpixels
    ind = [];
    return;
end

end

oldsegs = segs{ind(1)};

for k = 1:numel(oldsegs)
    if (numel(oldsegs{k})==numel(ind)) && oldsegs{k}(1)==ind(1) && all(oldsegs{k}==ind)
        ind = [];
        return;
    end
end

segs{ind(1)}{end+1} = ind;

         


