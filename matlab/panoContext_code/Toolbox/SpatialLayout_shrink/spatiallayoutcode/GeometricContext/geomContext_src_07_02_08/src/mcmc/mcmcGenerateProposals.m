function [smapNew, newseg, origseg, segadjmat, adjmat] =  ...
    mcmcGenerateProposals(pE, adjlist, smap, varargin)


% randomly split each segment further into connected components
%adjmat2 = curradjmat;
nsp = numel(smap);


repeatcount = 0;
while 1

    adjmat2 = zeros(nsp, nsp);
    for k = 1:size(adjlist, 1)    
        s1 = adjlist(k, 1);
        s2 = adjlist(k, 2); 
        if (smap(s1)==smap(s2)) && (rand(1) < pE(k))
            adjmat2(s1, s2) = 1;
            adjmat2(s2, s1) = 1;
        end    
    end

    % get segments by connected components
    smap2 = graphComponents(adjmat2);

    ncomponents = max(smap2);

    % select a new segment from the pool of fragments
    rcomp = ceil(rand(1)*ncomponents);

    % get indices of new segment
    s2ind = find(smap2==rcomp);
    
    % get index number of segment originally containing new segment
    origseg = smap(s2ind(1));

    % assign segment number to new segment
    smapNew = smap;
    if sum(smap==origseg) == numel(s2ind)    
        newseg = origseg;        
    else
        newseg = max(smapNew)+1;        
    end
    smapNew(s2ind) = newseg;
    
    % get sp adjacency matrix and segment adjacency matrix
    segadjmat = zeros(max(smapNew));
    adjmat = zeros(nsp, nsp); %curradjmat;
    for k = 1:size(adjlist, 1)
        s1 = adjlist(k, 1);
        s2 = adjlist(k, 2);
        seg1 = smapNew(s1);
        seg2 = smapNew(s2);
        if seg1==seg2
            adjmat(s1, s2) = 1;
            adjmat(s2, s1) = 1;
        else        
            segadjmat(seg1, seg2) = 1;
            segadjmat(seg2, seg1) = 1;
        end
    end    
    
    % make sure that reassigning  rcomp will not split the segment that it
    % is from 
    ind = find(smap == origseg);
    smapTmp = graphComponents(adjmat(ind, ind));    
    if max(smapTmp)<=2
        break;
    else
        repeatcount = repeatcount + 1;
        %disp(num2str([max(smapTmp) max(smapNew)]))
    end

%     if mod(repeatcount, 10)==0
%         disp(['repeats: ' num2str(repeatcount)])
%         if numel(varargin)==2
%             im = varargin{1};
%             segimage= varargin{2};
%             figure(1), displaySegmentGraph(im, smapNew(segimage), segadjmat);
%             figure(2), hist(smapNew, [1:max(smapNew)])
%             figure(3), hist(smapTmp, [1:max(smapTmp)])
%             figure(4), imagesc(label2rgb(smapNew(segimage))), axis image
%             figure(5), imagesc(label2rgb(smapTmp(segimage))), axis image
%             pause(1)
%         end
%     end
    
end            
