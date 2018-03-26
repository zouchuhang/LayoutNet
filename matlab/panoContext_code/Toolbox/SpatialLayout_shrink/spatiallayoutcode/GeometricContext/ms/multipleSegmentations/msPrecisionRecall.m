function pr = msPrecisionRecall(imsegs, labels, pg, cnum)

if ~exist('cnum', 'var') || isempty(cnum)
    cnum = 1;
end

conf = (0:0.01:1.0);
pcount = zeros(size(conf));
ncount = zeros(size(conf));

for k = 1:numel(imsegs)
    
    nind = labels{k}(:)==-1;
    pind = labels{k}(:)==1;
    
    tmppg = pg{k}(:, cnum);
    confim = tmppg(imsegs(k).segimage);
    pconf = confim(pind);
    nconf = confim(nind);
    
    npix = numel(confim);
    
    for c = 1:numel(conf)
        pconf = pconf(pconf>=conf(c));
        pcount(c) = pcount(c) + numel(pconf)/npix;
        nconf = nconf(nconf>=conf(c));
        ncount(c) = ncount(c) + numel(nconf)/npix;    
    end
        
end

pr.conf = conf;
pr.p = pcount ./ (pcount + ncount);
pr.r = pcount ./ (pcount(1));
pr.pcount = pcount;
pr.ncount = ncount;
pr.num = pcount(1);