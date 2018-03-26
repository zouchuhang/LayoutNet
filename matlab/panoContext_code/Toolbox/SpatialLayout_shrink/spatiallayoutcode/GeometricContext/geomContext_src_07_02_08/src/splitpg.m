function [pv, ph] = splitpg(pg)
% [pv, ph] = splitpg(pg)

if iscell(pg)
    for f = 1:numel(pg)
        pv{f} = [pg{f}(:, 1) sum(pg{f}(:, 2:6), 2) pg{f}(:, 7)];
        ph{f} = pg{f}(:, 2:6) ./ repmat(max(pv{f}(:, 2), 1E-6), [1 5]);
    end
else
    pv = [pg(:, 1) sum(pg(:, 2:6), 2) pg(:, 7)];
    ph = pg(:, 2:6) ./ repmat(max(pv(:, 2), 1E-6), [1 5]);    
end