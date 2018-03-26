function [vacc, hacc, vcm, hcm, pg] = ...
    testMultipleSegmentationsCV3(imsegs, labdata, segdata, maps, ...
    vclassifier, hclassifier, svclassifier, shclassifier, pvSP, phSP, ncv)
% [vacc, hacc, vcm, hcm] = testMultipleSegmentationsCV2(imsegs, labdata, segdata, maps, vclassifier, hclassifier, sclassifier, ncv)

pg = cell(numel(imsegs), 1);

for f = 1:numel(imsegs)
    
    if mod(f, 25)==0
        disp([num2str(f)])
    end
    
    nsp = imsegs(f).nseg;
    
    c = ceil(f/numel(imsegs)*ncv);

    pv{f} = pvSP{f};
    ph{f} = phSP{f};
    
    segs = cell(nsp, 1);
    
    for k = 1:size(labdata, 2)
        
        for s = 1:size(labdata{f, k}, 1)

            [segs, ind] = checksegs(segs, maps{f}(:, k), s);            
            
            if ~isempty(ind)
            
                vconf = test_boosted_dt_mc(vclassifier(c), labdata{f, k}(s, :));
                vconf = 1 ./ (1+exp(-vconf));
                vconf = vconf / sum(vconf); 

                hconf = test_boosted_dt_mc(hclassifier(c), labdata{f, k}(s, :));
                hconf = 1 ./ (1+exp(-hconf));
                hconf = hconf / sum(hconf);            

                svconf = test_boosted_dt_mc(svclassifier(c), segdata{f, k}(s, :));
                svconf = 1 ./ (1+exp(-svconf));           
                
                shconf = test_boosted_dt_mc(shclassifier(c), segdata{f, k}(s, :));
                shconf = 1 ./ (1+exp(-shconf));                 

                pv{f}(ind, :) = pv{f}(ind, :) + repmat(vconf*svconf, numel(ind), 1);
                ph{f}(ind, :) = ph{f}(ind, :) + repmat(hconf*shconf, numel(ind), 1);
                
            end
            
        end
                
    end
        
    pv{f} = pv{f} ./ max(repmat(sum(pv{f}, 2), 1, size(pv{f}, 2)), 1E-9);
    ph{f} = ph{f} ./ max(repmat(sum(ph{f}, 2), 1, size(ph{f}, 2)), 1E-9);    
    pg{f} = [pv{f}(:, 1)  repmat(pv{f}(:, 2), 1, 5).*ph{f}  pv{f}(:, 3)];      
         
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

         


