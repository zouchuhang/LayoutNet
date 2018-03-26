function [super, sub] = supersets(map1, map2)
% Return the segments (super) in map1 that are supersets of segments (sub)
% in map2. 
%
% Input:
%  map1(na) - a mapping from atoms to segments, such that 
%             map1(i)==map1(j)<-->i and j in same segment
%  map2(na) - a mapping from atoms to segments, such that 
%             map2(i)==map2(j)<-->i and j in same segment 
%
% Output:
%   super(nsuper) - segment in map1 that is superset of sub in map2
%   sub{nsuper}   - array of segments in map2 that are subset of a segment
%                   in map 1

super = [];
sub = [];

nseg1 = max(map1);
nseg2 = max(map2);

count2 = zeros(nseg2, 1);
for k = 1:nseg2
    count2(k) = sum(map2==k);
end

nsuper  = 0;

for k = 1:nseg1
    ind = find(map1==k);
    map2segs = unique(map2(ind), []);
    
    if numel(map2segs)>1 
        issuperset = 1;
        for s2 = 1:numel(map2segs)
            tmpcnt = sum(map2(ind)==map2segs(s2));
            if tmpcnt~=count2(map2segs(s2))
                issuperset = 0;
                break;
            end
        end
        if issuperset
            nsuper = nsuper+1;
            super(nsuper) = k;
            sub{nsuper} = map2segs;
        end
    end
end
     