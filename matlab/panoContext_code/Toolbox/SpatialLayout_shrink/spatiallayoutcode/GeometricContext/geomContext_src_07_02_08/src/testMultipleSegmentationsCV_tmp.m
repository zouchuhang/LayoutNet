function [vacc, hacc, vcm, hcm, pg, pg2] = ...
    testMultipleSegmentationsCV_tmp(imsegs, labdata, segdata, maps, ...
    vclassifier, hclassifier, sclassifier, pvSP, phSP, ncv)
% [vacc, hacc, vcm, hcm] = testMultipleSegmentationsCV_tmp(imsegs, labdata,
%                segdata, maps, vclassifier, hclassifier, sclassifier, ncv)

pg = cell(numel(imsegs), 1);

for f = 1:numel(imsegs)
    
    if mod(f, 1)==0
        disp([num2str(f)])
    end
    
    nsp = imsegs(f).nseg;
    
    c = ceil(f/numel(imsegs)*ncv);
   
    pg3{f} = [pvSP{f}(:, 1) repmat(pvSP{f}(:, 2), 1, 5).*phSP{f} pvSP{f}(:, 3)];   
    pg{f} = zeros(nsp, 7);
    
    segs = cell(nsp, 1);

    %count = zeros(nsp, 1);

    tmppg = zeros([size(pg{f}) size(labdata, 2)+1]);
    tmppg(:, :, end) = [pvSP{f}(:, 1) repmat(pvSP{f}(:, 2), 1, 5).*phSP{f} pvSP{f}(:, 3)];
     
    for k = 1:size(labdata, 2)
        
        for s = 1:size(labdata{f, k}, 1)

            %[segs, ind] = checksegs(segs, maps{f}(:, k), s);            
            ind = find(maps{f}(:, k)==s);
            
            if ~isempty(ind) 
            
                %count(ind) = count(ind) + 1;
                
                vconf = test_boosted_dt_mc(vclassifier(c), labdata{f, k}(s, :));
                vconf = 1 ./ (1+exp(-vconf));
                vconf = vconf / sum(vconf); 

                hconf = test_boosted_dt_mc(hclassifier(c), labdata{f, k}(s, :));
                hconf = 1 ./ (1+exp(-hconf));
                hconf = hconf / sum(hconf);            

                sconf = test_boosted_dt_mc(sclassifier(c), segdata{f, k}(s, :));
                sconf = 1 ./ (1+exp(-sconf));           

                pgs = [vconf(1) vconf(2)*hconf vconf(3)]*sconf;
                pg3{f}(ind, :) = pg3{f}(ind, :) + repmat(pgs, numel(ind), 1);
                
                pgs = [vconf(1) vconf(2)*hconf vconf(3)];
                tmppg(ind, :, k) = repmat(pgs, numel(ind), 1);
                
            end
            
        end
                
    end
   
    
    pv = cat(2, tmppg(:, 1, :), sum(tmppg(:, 2:6, :), 2), tmppg(:, 7, :));
    [val, tmppg_max] = max(pv, [], 2);    
    ind = find(tmppg_max==2);
    [val, tmppg_maxh] = max(tmppg(:, 2:6, :), [], 2);
    tmppg_max(ind) = 1+tmppg_maxh(ind);
    tmppg_max = squeeze(tmppg_max);
    
    sconf = zeros(size(tmppg));
    sconf(:, :, end) = 1;
    
    for k = 1:size(labdata, 2)
        for s = 1:size(labdata{f, k}, 1)
            ind = find(maps{f}(:, k)==s);
            sconf(ind, :, k) = repmat(...
                mean(all(tmppg_max(ind, [1:k-1 k+1:end])== ...
                repmat(tmppg_max(ind(1), [1:k-1 k+1:end]), [numel(ind) 1]), 1)), ...
                [numel(ind) 7]);
            %disp(num2str([sconf(ind(1), 1, k) tmppg_max(ind(1), k) numel(ind)]))
        end
    end

    pg{f} = sum(tmppg .* sconf, 3);    
    pg{f} = pg{f} ./ max(repmat(sum(pg{f}, 2), 1, size(pg{f}, 2)), 0.00001);    
    
    pg2{f} = sum(tmppg, 3);
    pg2{f} = pg2{f} ./ max(repmat(sum(pg2{f}, 2), 1, size(pg{f}, 2)), 0.00001);  
    
    pg3{f} = pg3{f} ./ max(repmat(sum(pg3{f}, 2), 1, size(pg3{f}, 2)), 0.00001);  
    
    [vacc, hacc] = mcmcProcessResult(imsegs(f), pg(f));
    [vacc2, hacc2] = mcmcProcessResult(imsegs(f), pg2(f));
    [vacc3, hacc3] = mcmcProcessResult(imsegs(f), pg3(f));
    %disp(num2str([vacc3 hacc3]))
    disp(['v: ' num2str([vacc2 vacc3 vacc]) '   h: ' num2str([hacc2 hacc3 hacc])])
    
%    [vacc, hacc] = testMultipleSegmentationsCV2(imsegs(f), labdata(f, :), segdata(f, :), maps(f), ...
%            vclassifier(c), hclassifier(c), sclassifier(c), pvSP(f), phSP(f), 1)
    
end

[vacc, hacc, vcm, hcm] = mcmcProcessResult(imsegs, pg);
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [segs, ind] = checksegs(segs, map, s)
% Checks whether this segment has been seen before

ind = find(map==s);

if numel(ind)==1 % already accounted for by superpixels
    ind = [];
    return;
end

oldsegs = segs{ind(1)};

for k = 1:numel(oldsegs)
    if (numel(oldsegs{k})==numel(ind)) && all(oldsegs{k}==ind)
        ind = [];
        return;
    end
end

segs{ind(1)}{end+1} = ind;

         


