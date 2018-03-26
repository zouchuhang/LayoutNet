function [vacc, hacc, vcm, hcm] = ...
    testMultipleSegmentationsCV(imsegs, labdata, segdata, maps, vclassifier, hclassifier, sclassifier, ncv)

for f = 1:numel(imsegs)
    
    c = ceil(f/numel(imsegs)*ncv);

    pg{f} = zeros(imsegs(f).nseg, 7);    
    
    for k = 1:size(labdata, 2)
        
        for s = 1:size(labdata{f, k}, 1)

            vconf = test_boosted_dt_mc(vclassifier(c), labdata{f, k}(s, :));
            vconf = 1 ./ (1+exp(-vconf));
            vconf = vconf / sum(vconf); 
            
            hconf = test_boosted_dt_mc(hclassifier(c), labdata{f, k}(s, :));
            hconf = 1 ./ (1+exp(-hconf));
            hconf = hconf / sum(hconf);            
            
            sconf = test_boosted_dt_mc(sclassifier(c), segdata{f, k}(s, :));
            sconf = 1 ./ (1+exp(-sconf));           
            
            ind = find(maps{f}(:, k)==s);
            
            pg{f}(ind, :) = pg{f}(ind, :) + repmat([vconf(1) vconf(2)*hconf vconf(3)], numel(ind), 1)*sconf;
            
        end
                
    end
        
    pg{f} = pg{f} ./ max(repmat(sum(pg{f}, 2), 1, size(pg{f}, 2)), 0.00001);    
        
end

[vacc, hacc, vcm, hcm] = mcmcProcessResult(imsegs, pg);
        
            
         


