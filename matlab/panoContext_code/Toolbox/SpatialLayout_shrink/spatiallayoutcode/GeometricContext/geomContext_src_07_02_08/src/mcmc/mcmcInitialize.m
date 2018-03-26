function [pvSP, phSP, pE, smap, adjmat] = ...
    mcmcInitialize(spdata, edata, adjlist, imsegs, vclassifierSP, hclassifierSP, eclassifier, ecal, inittype)

% probability of superpixel main labels
pvSP = test_boosted_dt_mc(vclassifierSP, spdata);
pvSP = 1 ./ (1+exp(-pvSP));
pvSP = pvSP ./ repmat(sum(pvSP, 2), 1, size(pvSP, 2));
[tmp, vmax] = max(pvSP, [], 2);

% probability of superpixel sub labels
phSP = test_boosted_dt_mc(hclassifierSP, spdata);
phSP = 1 ./ (1+exp(-phSP));
phSP = phSP ./ repmat(sum(phSP, 2), 1, size(phSP, 2));
[tmp, hmax] = max(phSP, [], 2);

% edge probability
pE = test_boosted_dt_mc(eclassifier, edata);
pE = 1 ./ (1+exp(ecal(1)*pE+ecal(2)));

%disp(num2str([min(pE) max(pE) mean(pE)]))

if exist('inittype') && strcmp(lower(inittype), 'none')
    smap = [];
    adjmat = [];
    return;
end


adjmat = zeros(size(imsegs.adjmat));
% random initialization according to pE
if exist('inittype') && strcmp(lower(inittype), 'random')
    disp('random')
    for k = 1:size(adjlist, 1)    
        s1 = adjlist(k, 1);
        s2 = adjlist(k, 2); 
        if rand(1) < pE(k)
            adjmat(s1, s2) = 1;
            adjmat(s2, s1) = 1;
        end    
    end
% create graph based on most likely edges
elseif exist('inittype') && strcmp(lower(inittype), 'max')
    for k = 1:size(adjlist, 1)    
        s1 = adjlist(k, 1);
        s2 = adjlist(k, 2); 
        if pE(k) > 0.99
            adjmat(s1, s2) = 1;
            adjmat(s2, s1) = 1;
        end    
    end
% create initial graph based on most likely superpixel labels
else 
    for k = 1:size(adjlist,1)
        s1 = adjlist(k, 1);
        s2 = adjlist(k, 2);
        if (vmax(s1)==vmax(s2)) && (vmax(s1)~=2 || (hmax(s1)==hmax(s2)))
            adjmat(s1, s2) = 1;
            adjmat(s2, s1) = 1;
        end
    end
end



% get segments by connected components
smap = graphComponents(adjmat);

%disp(num2str(max(smap)))

% lmap = zeros(max(smap), 1);
% for k = 1:numel(lmap)
%     ind = find(smap==k);
%     lmap(k) = vmax(ind(1)) + (vmax(ind(1))==2)*hmax(ind(1))/10;
% end
