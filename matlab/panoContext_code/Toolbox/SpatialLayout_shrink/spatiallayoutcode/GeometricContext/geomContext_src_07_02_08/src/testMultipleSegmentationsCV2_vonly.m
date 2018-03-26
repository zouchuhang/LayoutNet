function [vacc, vcm, pg] = ...
    testMultipleSegmentationsCV2(imsegs, labdata, segdata, maps, ...
    vclassifier, sclassifier, pvSP, ncv)
% [vacc, hacc, vcm, hcm] = testMultipleSegmentationsCV2(imsegs, labdata, segdata, maps, vclassifier, hclassifier, sclassifier, ncv)

pg = cell(numel(imsegs), 1);
allcount = 0;

for f = 1:numel(imsegs)
    
    if mod(f, 50)==0
        disp([num2str(f) ' ' num2str(allcount/(f-1))])
    end
    
    nsp = imsegs(f).nseg;
    
    c = ceil(f/numel(imsegs)*ncv);

    pg{f} = [pvSP{f}];%(:, 1) repmat(pvSP{f}(:, 2), 1, 5).*phSP{f} pvSP{f}(:, 3)];   
    %pg{f} = zeros(nsp, 7);
    
    segs = cell(nsp, 1);

    count = zeros(nsp, 1);
    
    for k = 1:size(labdata, 2)
        
        for s = 1:size(labdata{f, k}, 1)

            [segs, ind] = checksegs(segs, maps{f}(:, k), s);            
            
            if ~isempty(ind)
            
                count(ind) = count(ind) + 1;
                
                vconf = test_boosted_dt_mc(vclassifier(c), labdata{f, k}(s, :));
                vconf = 1 ./ (1+exp(-vconf));
                vconf = vconf / sum(vconf); 
% 
%                 hconf = test_boosted_dt_mc(hclassifier(c), labdata{f, k}(s, :));
%                 hconf = 1 ./ (1+exp(-hconf));
%                 hconf = hconf / sum(hconf);            

                sconf = test_boosted_dt_mc(sclassifier(c), segdata{f, k}(s, :));
                sconf = 1 ./ (1+exp(-sconf));           

                %pgs = [vconf(1) vconf(2)*hconf vconf(3)]*sconf;
                pgs = vconf*sconf;

                pg{f}(ind, :) = pg{f}(ind, :) + repmat(pgs, numel(ind), 1);
            end
            
        end
                
    end
        
    pg{f} = pg{f} ./ max(repmat(sum(pg{f}, 2), 1, size(pg{f}, 2)), 0.00001);    
        
    allcount = allcount + mean(count);

    % convert to seven class format (which will be ignored)
    pg{f} = [pg{f}(:, 1) repmat(pg{f}(:, 2), [1 5])/5 pg{f}(:, 3)];
    
end

[vacc, hacc, vcm, hcm] = mcmcProcessResult(imsegs, pg);
        
allcount = allcount / numel(imsegs);


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

         


