function [maps, rscores, Z] = APPsp2regions(segDensity, spdata, nSegments)
% [maps, rscores, Z] = APPsp2regions(segDensity, spdata, nSegments)
% Clusters the superpixels into segments according to pairwise likelihoods
% given by segDensity
%
% Input:
%   segDensity: pairwise sp likelihoods structure
%   spdata: superpixel data and features structure
%   nSegments: array containing the number of segments for each
%              segmentation
% Ouput:
%   maps(nSp, nMaps): randomly generated clusterings according
%                     to the distribution given by densities and features
%   rscores: the average log likelihood for each region
%   Z: the pairwise superpixel likelihood matrix
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

features = spdata.features;

nSp = size(features, 1);
nFeatures = size(features, 2);
nDist = nSp*(nSp-1)/2;
nMaps = length(nSegments);

Z = APPgetPairwiseSuperpixelLikelihoods(features, segDensity);

% convert Z to P(y1=y2 | |x1-x2|)
Z = log(1./(1+exp(-Z)));

maps = zeros(nSp, nMaps);
for m = 1:nMaps
    
    num_clusters = nSegments(m);
    if num_clusters > nSp
        num_clusters = nSp;
    end
    % get random ordering of indices
    pind = randperm(nSp);
       
    % assign cluster ids to first num_clusters clusters    
    for i = 1:num_clusters
        maps(pind(i), m) = i;       
    end

    % assign remaining clusters
    % go forwards ...
    for i = num_clusters+1:length(pind)
        % compute P(c(i) = k | x) for each k
        f = zeros(num_clusters, 1);
        likelihoods = Z(pind(1:i-1), pind(i));            
        for k = 1:num_clusters
            sameinds = find(maps(pind(1:i-1), m)==k);         
            f(k) = mean(likelihoods(sameinds));
        end 
       
        [tmp, k] = max(f);
        
        % assign cluster and add log(P(ci = k | x, c))
        maps(pind(i), m) = k;
        %fact = repmat(-1, [i-1 1]);
        %fact(find(maps(pind(1:i-1), m)==k))=1;

    end
    
    %... and go backwards now using all cluster assignments
    for pp = 2:2
        %for i = num_clusters+1:length(pind)
        for i = 1:length(pind) % allow initial assignments to be reassigned
            % compute P(c(i) = k | x) for each k
            f = zeros(num_clusters, 1);
            likelihoods = Z(:, pind(i));            
            for k = 1:num_clusters
                sameinds = find(maps(:, m)==k);          
                if length(sameinds) > 0
                    f(k) = mean(likelihoods(sameinds));
                end
            end 
	     
            [tmp, k] = max(f);
            
            % assign cluster and add log(P(ci = k | x, c))
            maps(pind(i), m) = k;
            %fact = repmat(-1, [i-1 1]);
            %fact(find(maps(pind(1:i-1), m)==k))=1;        
        end
    end
    
    currmap = maps(:, m);
    nc = max(currmap);
    for k = 1:nc
        ninds = length(find(currmap==k));
        while ninds == 0 & k < nc
            inds = find(currmap>k);
            currmap(inds) = currmap(inds)-1;
            ninds = length(find(currmap==k));
            nc = max(currmap);
        end
    end
    maps(:, m) = currmap;
    
    if nargout > 1
        rscores{m} = zeros(nc, 1);
        for k = 1:nc
            ind = find(currmap==k);
            npairs = length(ind)*(length(ind)-1);
            % average of average log-likelihoods
            rscores{m}(k) = 1./(1+exp(-sum(sum(Z(ind, ind)))/npairs));
        end
    end

end

