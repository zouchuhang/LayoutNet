function [cimages, cnames] = APPclassifierOutput2confidenceImages(imsegs, conf_maps)
% Computes confidence maps from results of APPtestImage
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

nvclasses = length(conf_maps(1).vnames);
nhclasses = length(conf_maps(1).hnames);
nclasses = nvclasses+nhclasses;
for v = 1:nvclasses
    cnames{v} = ['v' conf_maps(1).vnames{v}];
end
for h = 1:nhclasses
    cnames{nvclasses+h} = ['h' conf_maps(1).hnames{h}];
end


cimages = cell(length(imsegs), 1);

for f = 1:length(imsegs)
 
    for n = 1:nclasses
        cim{n} = single(ones(imsegs(f).imsize)/nclasses);
    end
        
%     for s = 1:imsegs(f).nseg
%         
%         ind = find(imsegs(f).segimage==s);
%         
%         for v = 1:nvclasses
%             cim{v}(ind) = single(conf_maps(f).vmap(s, v));            
%         end
%         for h = 1:nhclasses
%             cim{nvclasses+h}(ind) = single(conf_maps(f).hmap(s, h));
%         end
%         
%     end
    
    cimages{f} = single(zeros([imsegs(f).imsize nclasses]));
    
    vind = find(strcmp(cnames, 'v090'));      
    
    for n = 1:nclasses
        if n <= nvclasses
            %cimages{f}(:, :, n) = cim{n};
            tmpmap = conf_maps(f).vmap(:, n);
            cimages{f}(:, :, n) = tmpmap(imsegs(f).segimage);
        else
            % multiply P(vclass|X) by P(vertical|X)            
            tmpmap = conf_maps(f).hmap(:, n-nvclasses);
            %cimages{f}(:, :, n) = cim{n}.*cim{vind}; 
            cimages{f}(:, :, n) = tmpmap(imsegs(f).segimage).*cimages{f}(:, :, vind);
        end
    end    
end
    