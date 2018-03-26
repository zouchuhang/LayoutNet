function [prv, prh] = pg2prcurve(pg, imsegs)
% For varying levels of confidence, compute roc curves (tp vs fp) for the
% vertical labels and the horizontal labels
% labels(num_im).{vert_labels(h, w), vert_conf(h, w), horz_labels(h, w),
% horz_conf(h, w)

[pv, ph] = splitpg(pg);

[y, c, w] = initData(pv, {imsegs(:).vert_labels}, {imsegs(:).npixels});
prv = computePR(y, c, w);

[y, c, w] = initData(ph, {imsegs(:).horz_labels}, {imsegs(:).npixels});
prh = computePR(y, c, w);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pr = computePR(y, c, w)

for v = 1:numel(y)
    [c{v}, ind] = sort(c{v}, 'descend');
    y{v} = y{v}(ind);
    w{v} = w{v}(ind);
    w{v} = w{v} / sum(w{v});
    pr(v).conf = c{v};
    pr(v).p = cumsum((y{v}==1).*w{v}) ./ cumsum(w{v});
    pr(v).r = cumsum((y{v}==1).*w{v}) / sum((y{v}==1).*w{v});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, c, w] = initData(p, lab, npix)

ndata = 0;
for f = 1:numel(p)
    ndata = ndata + sum(lab{f}>0);
    npix{f} = npix{f} / sum(npix{f});
end

nclasses = size(p{1}, 2);

for v = 1:nclasses
    y{v} = zeros(ndata, 1);
    c{v} = zeros(ndata, 1);
    w{v} = zeros(ndata, 1);

    n = 0;
    for f = 1:numel(p)
        ind = find(lab{f}>0);
        nsp = numel(ind);
        y{v}(n+1:n+nsp) = (lab{f}(ind)==v);
        c{v}(n+1:n+nsp) = p{f}(ind, v);
        w{v}(n+1:n+nsp) = npix{f}(ind);
        n = n + nsp;
    end
end

