function [boundmap, perim] = mcmcGetSuperpixelBoundaries_fast(imsegs)
% boundmap{nimages}{nseg, nseg}(npixidx) - boundaries between pairs of segs
% perim{nimages}(nseg, nseg) - number of pixels in boundary
% nimages == 1, boundmap{nseg, nseg}(npixidx), perim(nseg, nseg)

for f = 1:numel(imsegs)
    
    nseg = imsegs(f).nseg;	
    segimage = imsegs(f).segimage;
    [imh, imw] = size(segimage);
    
    % get adjacency
    dx = segimage ~= segimage(:,[2:end end]);
    dy = segimage ~= segimage([2:end end], :);    
    
    ind1 = find(dy);
    ind2 = ind1 + 1;
    s1 = segimage(ind1);
    s2 = segimage(ind2);
    ind3 = find(dx);
    ind4 = ind3 + imh;
    s3 = segimage(ind3);
    s4 = segimage(ind4);
    
    % get boundaries
    ind = [ind1 ; ind3];
    s12 = [[min([s1 s2], [], 2) max([s1 s2], [], 2)] ; ...
        [min([s3 s4], [], 2) max([s3 s4], [], 2)]];
    perim{f} = zeros(nseg, nseg, 'uint16');
    s1 = s12(:, 1);  s2 = s12(:, 2);
    for k = 1:numel(ind)
        perim{f}(s1(k), s2(k)) = perim{f}(s1(k), s2(k))+1;
    end
    bndind = find(perim{f}>0);
    boundmap{f}= cell(nseg, nseg);
    for k = bndind'
        ts1 = mod(k-1, nseg)+1;
        ts2 = floor((k-1)/nseg)+1;
        boundmap{f}{k} = ind(s1==ts1 & s2==ts2);
    end    
       
end
    
if f==1
    boundmap = boundmap{1};
    perim = perim{1};
end
