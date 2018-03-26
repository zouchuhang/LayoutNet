function pg = msTestRf(imsegs, data, maps, dtlab, dtseg, normalize)
% [vacc, hacc, vcm, hcm] = testMultipleSegmentationsCV2(imsegs, labdata, segdata, maps, vclassifier, hclassifier, sclassifier, ncv)

if ~exist('normalize', 'var') || isempty(normalize)
    normalize = 1;
end

pg = cell(numel(imsegs), 1);

nclasses = size(dtlab(1).classcount,2);

for f = 1:numel(imsegs)
            
    disp(num2str(f))
    
    nsp = imsegs(f).nseg;   

    pg{f} = zeros(nsp, nclasses);
    
    
    segs = cell(nsp, 1);

    if ~normalize
        sum_s = zeros(nsp, 1);
    end
    
    for k = 1:size(data, 2)
        
        isvalid = false(size(data{f, k}, 1), 1);
        for s = 1:size(data{f, k}, 1)
            [segs, ind] = checksegs(segs, maps{f}(:, k), s);            
            isvalid(s) = ~isempty(ind);
        end
            
        yconf_all = rfTreeVal(dtlab, data{f, k}(isvalid, :)); 
        sconf_all = rfTreeVal(dtseg, data{f, k}(isvalid, :));  
                  
        isvalid = find(isvalid);
        for s = 1:numel(isvalid)     
            
            ind = find(maps{f}(:, k)==isvalid(s));


            yconf = yconf_all(s, :);
            %yconf = 1 ./ (1+exp(-yconf));
            
            if normalize
                yconf = yconf / sum(yconf);        
            end

            sconf = sconf_all(s, 1);
            %sconf = 1 ./ (1+exp(-sconf));           
          
            pgs = yconf*sconf;
            
            if ~normalize
                sum_s(ind) = sum_s(ind) + sconf;
            end
            pg{f}(ind, :) = pg{f}(ind, :) + repmat(pgs, numel(ind), 1);
            
        end
                
    end
     
    if normalize
        pg{f} = pg{f} ./ max(repmat(sum(pg{f}, 2), 1, size(pg{f}, 2)), 0.00001);    
    else
        pg{f} = pg{f} ./ repmat(sum_s, [1 size(pg{f}, 2)]);
    end 
        
    
end
        

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

         


