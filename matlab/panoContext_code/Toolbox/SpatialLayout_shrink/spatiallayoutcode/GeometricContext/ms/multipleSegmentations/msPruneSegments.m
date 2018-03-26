%% Remove duplicate segments
function smaps2 = msPruneSegments(smaps)

smaps2 = smaps;
segment_list = {};

for m = 1:size(smaps, 2)
    for k = 1:max(smaps(:, m))
        ind = smaps(:, m)==k;
        [segment_list, isnew] = prune_add2list(segment_list, find(ind));
        if ~isnew
            smaps2(ind, m) = 0;
        end
    end
    uid = unique(smaps2(:, m));
    uid = setdiff(uid, 0);
    for k = 1:numel(uid)
        ind = smaps2(:, m)==uid(k);
        smaps2(ind, m) = k;
    end
end


function [seglist, isnew] = prune_add2list(seglist, ids)

isnew = true;
for k = 1:numel(seglist)
    seg = seglist{k};
    if numel(seg)==numel(ids) && all(seg==ids)
        isnew = false;
        return;
    end
end
seglist{end+1} = ids;