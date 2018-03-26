function pg = confidenceImages2pg(idx, cim)
% pg = confidenceImages2pg(idx, cim)
% Computes average likelihood for regions given by pixel indices
%

NO_CELL = 0;

if ~iscell(idx{1})
    idx = {idx};
    NO_CELL = 1;    
end
if ~iscell(cim)
    cim = {cim};    
    NO_CELL = 1;     
end

for f = 1:numel(idx)    

    pg{f} = zeros(numel(idx{f}), size(cim{f}, 3));    
    
    for c = 1:size(cim{f}, 3)
        
        tmpc = cim{f}(:, :, c);
        
        for k = 1:numel(idx{f})
            pg{f}(k, c) = sum(tmpc(idx{f}{k})) / numel(idx{f}{k});
        end
    end
end
    
if NO_CELL
    pg = pg{f};
end
    