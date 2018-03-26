function [vacc, hacc, vcm, hcm, pg] = ...
    testSingleSegmentationsCV(imsegs, labdata, maps, vclassifier, hclassifier, ncv)
% [vacc, hacc, vcm, hcm, pg] = testSingleSegmentationsCV(
%                       imsegs, labdata, maps, vclassifier, hclassifier, ncv)

pg = cell(numel(imsegs), 1);

for f = 1:numel(imsegs)
       
    nsp = imsegs(f).nseg;
    
    c = ceil(f/numel(imsegs)*ncv);

    pg{f} = zeros(nsp, 7);  
    
    segs = cell(nsp, 1);

    for s = 1:size(labdata{f}, 1)

        ind = find(maps{f}==s);
        
        if ~isempty(ind)
                           
            vconf = test_boosted_dt_mc(vclassifier(c), labdata{f}(s, :));
            vconf = 1 ./ (1+exp(-vconf));
            vconf = vconf / sum(vconf); 

            hconf = test_boosted_dt_mc(hclassifier(c), labdata{f}(s, :));
            hconf = 1 ./ (1+exp(-hconf));
            hconf = hconf / sum(hconf);                     

            pgs = [vconf(1) vconf(2)*hconf vconf(3)];

            pg{f}(ind, :) = repmat(pgs, numel(ind), 1);
            
        end
                
    end                 
    
end

[vacc, hacc, vcm, hcm] = mcmcProcessResult(imsegs, pg);
        



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [segs, ind] = checksegs(segs, map, s)
% Checks whether this segment has been seen before

ind = find(map==s);

oldsegs = segs{ind(1)};

for k = 1:numel(oldsegs)
    if (numel(oldsegs{k})==numel(ind)) && all(oldsegs{k}==ind)
        ind = [];
        return;
    end
end

segs{ind(1)}{end+1} = ind;

         


