function splabels = msLabelMap2Sp(lmap, smap)

nocell = false;
if ~iscell(lmap)
    nocell = true;
    lmap = {lmap};
    smap = {smap};
end

% if ~exist('maxlab', 'var') || isempty(maxlab)
%     maxlab = 0;
%     for f = 1:numel(lmap)
%         maxlab = max(maxlab, max(lmap{f}(:)));
%     end
% end

for f = 1:numel(lmap)

    nseg = max(smap{f}(:));
    stats = regionprops(smap{f}, 'PixelIdxList');
    idx = {stats.PixelIdxList};
   
    splabels{f} = zeros(nseg, 1);
    for s = 1:nseg
        splabels{f}(s) = mode(double(lmap{f}(idx{s})));
    end
end
   
if nocell
    splabels = splabels{1};
end