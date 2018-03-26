function cinds = APPgetSpIndsOld(map, imsegs);
% Get pixel indices for each superpixel

nc = max(map(:));
cinds(nc) = struct('inds', [], 'count', 0);

t1 = cputime;
for c = 1:nc
    spind = find(map==c);
    nind = sum(imsegs.npixels(spind));
    cinds(c).inds = zeros(nind, 1);
    cinds(c).count = 0;
end

t2 = cputime;
nseg = imsegs.nseg;
for k = 1:nseg
    ind = find(imsegs.segimage==k);
    c = map(k);
    cinds(c).inds(cinds(c).count+1:cinds(c).count+length(ind)) = ind;
    cinds(c).count = cinds(c).count + length(ind);
end

for c = 1:length(cinds)
    cinds(c).inds = cinds(c).inds(1:cinds(c).count);
end
%segimage = imsegs.segimage;
%for k = 1:numel(segimage)
%    c = map(segimage(k));
%    cinds(c).count = cinds(c).count + 1;
%    cinds(c).inds(cinds(c).count) = k;
%end

