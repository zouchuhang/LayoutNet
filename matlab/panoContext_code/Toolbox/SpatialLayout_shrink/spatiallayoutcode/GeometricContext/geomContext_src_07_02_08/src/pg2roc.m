function [vroc, hroc, vcm_s, hcm_s, vcm_p, hcm_p] = ...
    pg2roc(labels, confmaps, gtsegs)
% For varying levels of confidence, compute roc curves (tp vs fp) for the
% vertical labels and the horizontal labels
% labels(num_im).{vert_labels(h, w), vert_conf(h, w), horz_labels(h, w),
% horz_conf(h, w)


nvclasses = numel(gtsegs(1).vert_names);
nhclasses = numel(gtsegs(1).horz_names);

confidences = (0.0:0.02:1.0)';

vroc.conf = confidences;
vroc.fp = zeros(length(confidences), 1);
vroc.tp = zeros(length(confidences), 1);
vroc.count = 0;

hroc.conf = confidences;
hroc.fp = zeros(length(confidences), 1);
hroc.tp = zeros(length(confidences), 1);
hroc.count = 0;

for i =1:length(gtsegs)
    if ~strcmp(labels(i).imname, gtsegs(i).imname)
        disp(['name mismatch: ' num2str(i) ' ' gtsegs(i).imname ' ' labels(i).imname]);
    end
end

for i = 1:length(gtsegs)  
    for n = 1:numel(gtsegs(i).vert_labels)
        if gtsegs(i).vert_labels(n)~=0
            vroc.count = vroc.count + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
            if gtsegs(i).horz_labels(n)~=0
                hroc.count = hroc.count + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
            end
        end
    end
end

for c = 1:length(confidences)
    conf = confidences(c);    
    
    for i = 1:length(gtsegs)

        %vlab = repmat({'---'}, gtsegs(i).nseg, 1);
        %hlab = repmat({'---'}, gtsegs(i).nseg, 1);        
        
        %ind = find(gtsegs(i).vert_labels~=0);
        %vlab(ind) = gtsegs(i).vert_names(gtsegs(i).vert_labels(ind));            
        %ind = ind(find(gtsegs(i).horz_labels(ind)~=0));
        %hlab(ind) = gtsegs(i).horz_names(gtsegs(i).horz_labels(ind)); 
        vlab = gtsegs(i).vert_labels;
        hlab = gtsegs(i).horz_labels;
        
        
        for n = 1:numel(vlab)
            
            if (vlab(n)~=0) & (max(confmaps(i).vmap(n, :)) >= conf)                
                if (confmaps(i).vmap(n, vlab(n)) >= conf)                        
                    vroc.tp(c) = vroc.tp(c) + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
                end
                if max(confmaps(i).vmap(n, [1:vlab(n)-1 vlab(n)+1:end])) >= conf
                    vroc.fp(c) = vroc.fp(c) + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
                end
            end
            
            if (hlab(n)~=0) & (max(confmaps(i).hmap(n, :)) >= conf)                
                if (confmaps(i).hmap(n, hlab(n)) >= conf)                        
                    hroc.tp(c) = hroc.tp(c) + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
                end
                if max(confmaps(i).hmap(n, [1:hlab(n)-1 hlab(n)+1:end])) >= conf
                    hroc.fp(c) = hroc.fp(c) + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
                end
            end
            
        end
    end
end

vroc.tp = vroc.tp / vroc.count;
vroc.fp = vroc.fp / vroc.count;
hroc.tp = hroc.tp / hroc.count;
hroc.fp = hroc.fp / hroc.count;

% get confusion matricies
vnames = gtsegs(i).vert_names;
hnames = gtsegs(i).horz_names;

vcm_s.names = vnames;
hcm_s.names = hnames;
vcm_p.names = vnames;
hcm_p.names = hnames;
nvn = length(vnames); 
nhn = length(hnames);


vmat_s = zeros(nvn, nvn);
hmat_s = zeros(nhn, nhn);
vmat_p = zeros(nvn, nvn);
hmat_p = zeros(nhn, nhn);

for i = 1:length(gtsegs)

    vci = gtsegs(i).vert_labels;
    hci = gtsegs(i).horz_labels;
   
    
    for n = 1:numel(vci)

        if vci(n)~=0           
            v1 = vci(n);
            v2 = find(strcmp(labels(i).vert_labels(n), vnames));
            vmat_s(v1, v2) = vmat_s(v1, v2)+1;
            vmat_p(v1, v2) = vmat_p(v1, v2) + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
        end
                          
        if hci(n)~=0           
            h1 = hci(n);
            h2 = find(strcmp(labels(i).horz_labels(n), hnames));
            hmat_s(h1, h2) = hmat_s(h1, h2)+1;
            hmat_p(h1, h2) = hmat_p(h1, h2) + gtsegs(i).npixels(n)/prod(gtsegs(i).imsize);
        end     

    end
end

vcm_s.acc = 0;
vcm_p.acc = 0;
for i = 1:nvn    
    vcm_s.count(i) = sum(vmat_s(i, :));
    vcm_s.asgncount(i) = sum(vmat_s(:, i));
    vcm_s.mat(i, :) = vmat_s(i, :) / vcm_s.count(i);        
    vcm_p.count(i) = sum(vmat_p(i, :));
    vcm_p.asgncount(i) = sum(vmat_p(:, i));
    vcm_p.mat(i, :) = vmat_p(i, :) / vcm_p.count(i);   
end
for i = 1:nvn
    vcm_s.acc = vcm_s.acc + vcm_s.count(i)/sum(vcm_s.count)*vcm_s.mat(i,i);
    vcm_p.acc = vcm_p.acc + vcm_p.count(i)/sum(vcm_p.count)*vcm_p.mat(i,i);
end

hcm_s.acc = 0;
hcm_p.acc = 0;
for i = 1:nhn    
    hcm_s.count(i) = sum(hmat_s(i, :));
    hcm_s.asgncount(i) = sum(hmat_s(:, i));
    hcm_s.mat(i, :) = hmat_s(i, :) / hcm_s.count(i);        
    hcm_p.count(i) = sum(hmat_p(i, :));
    hcm_p.asgncount(i) = sum(hmat_p(:, i));
    hcm_p.mat(i, :) = hmat_p(i, :) / hcm_p.count(i);            
end
for i = 1:nhn
	hcm_s.acc = hcm_s.acc + hcm_s.count(i)/sum(hcm_s.count)*hcm_s.mat(i,i);
    hcm_p.acc = hcm_p.acc + hcm_p.count(i)/sum(hcm_p.count)*hcm_p.mat(i,i);   
end