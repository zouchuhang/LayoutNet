function smaps =  msCreateMultipleSegmentations(pE, adjlist, nsp, nsegall)
% 1) Randomly select superpixel s1, then randomly selects different superpixel
% s2 within same segment (if one exists); remove s1,s2 from s
% 2) Then, for randomly ordered i:
%     if si is adjacent to s1, assign si to s1 with probability pE(si, s1)
%     if si is adjacent to s2, assign si to s2 with probability pE(si, s2)
%     remove si from s
% 3) Repeat (2) until s is empty
%
% Note: this differs from the published version of the algorithm in that
% duplicate segments are removed (number is replaced with 0 in smaps)

% randomly split each segment further into connected components
%adjmat2 = curradjmat;
nmaps = numel(nsegall);

smaps = zeros(nsp, nmaps);

adjmat = zeros(nsp, nsp);
for k = 1:size(adjlist, 1)
    s1 = adjlist(k, 1);
    s2 = adjlist(k, 2);
    adjmat(s1, s2) = k;
    adjmat(s2, s1) = k;
end
    
nadj = zeros(nsp, 1);
for k = 1:nsp
    adj{k} = find(adjmat(k, :));
    nadj(k) = numel(adj{k});
end
for k = 1:nsp
    normalization = 1;%2./(nadj(k)+nadj(adj{k}));
    pOffAll{k} = (1-pE(adjmat(k, adj{k}))).^normalization + 1E-10;
    pOnAll{k} = pE(adjmat(k, adj{k})).^normalization + 1E-10;     
end

for m = 1:nmaps
    rind = randperm(nsp);
    smap = zeros(nsp, 1);
    nseg = nsegall(m);                    
    
    nseg = min(nseg, nsp);
    smap(rind(1:nseg)) = (1:nseg);
    
    if nseg==nsp % more segments than superpixels: return on sp per segment
        smaps(:, m) = smap;
        continue;
    end
    
    rind(1:nseg) = [];   

    nleft = nsp;
    %tic
    while sum(nadj(rind))>0 % do until all possible sp are assigned
        for r = rind
            %p = -Inf*ones(nseg, 1);
            
            if all(smap(adj{r})>0)
                p = -Inf*ones(nseg, 1);
                for k = unique(smap(adj{r}))'
                    %if any(smap(adj{r})==k)                           
                        sOn = (smap(adj{r})==k);
                        sOff = ~sOn; %setdiff(find(smap(adj{r})>0), sOn);
                        pOn = prod(pOnAll{r}(sOn));
                        pOff = 1;
                        if ~isempty(sOff)
                            pOff = prod(pOffAll{r}(sOff));   
                        end
                        p(k) = log(pOn)+log(pOff);                            
                    %end                                        
                end
                [kval, kmax] = max(p);
                smap(r) = kmax; 
            else
                for k = unique(smap(adj{r}))' %1:nseg       
                    if k > 0
                        kadj = (smap(adj{r})==k);
                    %if sum(kadj)>0
                        pOn = prod(pOnAll{r}(kadj));
                        pOff = prod(pOffAll{r}(kadj));
                        if rand(1) < pOn / (pOn+pOff) + 0.05
                            smap(r) = k;
                        end
                    %end
                    end
                end
            end
        end

        rind(smap(rind)>0) = [];                
        
        % break if it seems stuck (can happen when pE is very low,
        % nsegments is small)
        if numel(nleft)>50 && all(nleft(end-49:end)==numel(rind))
            %disp(['break init early: ' num2str(numel(rind)) ' left'])
            %disp(sort(rind))
            break;
        end
        nleft(end+1) = numel(rind);
        
    end

    count = zeros(nseg, 1);
    for k = 1:nseg
        count(k) = sum(smap==k);
    end
    %toc
    %disp(['0: ' num2str(evaluateEdgeProb(adjlist, pE, smap))])
    %tic
    for t = 1:10
        lastmap = smap;
        rind = randperm(nsp);
        for r = rind
            if smap(r) > 0 && count(smap(r)) > 1
                p = -Inf*ones(nseg, 1);               

                for k = unique(smap(adj{r}))' %1:nseg
                    if k > 0
                    %if sum(smap(adj{r})==k)>0                           
                        sOn = smap(adj{r})==k;
                        sOff = ~sOn; %setdiff([1:nadj(r)], sOn);
                        pOn = sum(log(pOnAll{r}(sOn)));
                        pOff = sum(log(pOffAll{r}(sOff)));   
                        p(k) = pOn + pOff; 
                    end
                    %end
                end
                [kval, kmax] = max(p);
                if kmax~=smap(r)
                    count(smap(r)) = count(smap(r))-1;
                    smap(r) = kmax;
                end

            end
        end
        if all(lastmap==smap)
            break;
        end
    end 
    %disp(['1: ' num2str(evaluateEdgeProb(adjlist, pE, smap))])
    %toc
    %disp(num2str(t))
    smaps(:, m) = smap;
end

smaps = msPruneSegments(smaps);

%%

function p = evaluateEdgeProb(adjlist, pE, smap)
hasedge = zeros(size(adjlist, 1), 1);
for k = 1:size(adjlist, 1)
    if smap(adjlist(k, 1))==smap(adjlist(k, 2))
        hasedge(k)=1;
    end
end
p = sum(log(pE(hasedge==1)));
p = p+sum(log(1-pE(hasedge==0)));




   



