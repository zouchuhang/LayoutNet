function smaps =  generateMultipleSegmentationsGC(pE, adjlist, nsp, nsegall)
% smaps =  generateMultipleSegmentationsGC(pE, adjlist, nsp, nsegall)
% 
% Uses graph cuts to segment image, maximizing an approximation of log(pE)
% Different segments result from random superpixels initially being set to 
% segments.  

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
    pOffAll{k} = (1-pE(adjmat(k, adj{k}))).^normalization;
    pOnAll{k} = pE(adjmat(k, adj{k})).^normalization;     
end

for m = 1:nmaps
    rind = randperm(nsp);
    smap = zeros(nsp, 1);
    nseg = nsegall(m);
    
    nseg = min(nseg, nsp);
    smap(rind(1:nseg)) = (1:nseg);
    
    rind(1:nseg) = [];   

    nleft = nsp;
    
    while sum(nadj(rind))>0 % do until all possible sp are assigned
        for r = rind
            p = -Inf*ones(nseg, 1);
            
            if all(smap(adj{r})>0)
                p = repmat(-Inf, [nseg 1]);
                for k = 1:nseg
                    if any(smap(adj{r})==k)                           
                        sOn = find(smap(adj{r})==k);
                        sOff = setdiff(find(smap(adj{r})>0), sOn);
                        pOn = prod(pOnAll{r}(sOn));
                        pOff = 1;
                        if ~isempty(sOff)
                            pOff = prod(pOffAll{r}(sOff));   
                        end
                        p(k) = log(pOn)+log(pOff);                            
                    end                                        
                end
                [kval, kmax] = max(p);
                smap(r) = kmax; 
            else
                for k = 1:nseg       
                    kadj = find(smap(adj{r})==k);
                    if ~isempty(kadj)
                        pOn = prod(pOnAll{r}(kadj));
                        pOff = prod(pOffAll{r}(kadj));
                        if rand(1) < pOn / (pOn+pOff)
                            smap(r) = k;
                        end
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
    
    %disp(['0: ' num2str(evaluateEdgeProb(adjlist, pE, smap))])
    %tic
    for t = 1:10
        lastmap = smap;
        rind = randperm(nsp);
        for r = rind
            if smap(r) > 0
                p = repmat(-Inf, [nseg 1]);               
                if count(smap(r)) > 1
                    for k = 1:nseg
                        if any(smap(adj{r})==k)                           
                            sOn = find(smap(adj{r})==k);
                            sOff = setdiff([1:nadj(r)], sOn);
                            pOn = prod(pOnAll{r}(sOn));
                            pOff = prod(pOffAll{r}(sOff));   
                            p(k) = log(pOn)+log(pOff);                            
                        end
                    end
                    [kval, kmax] = max(p);
                    if kmax~=smap(r)
                        count(smap(r)) = count(smap(r))-1;
                        smap(r) = kmax;
                    end
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




function p = evaluateEdgeProb(adjlist, pE, smap)
hasedge = zeros(size(adjlist, 1), 1);
for k = 1:size(adjlist, 1)
    if smap(adjlist(k, 1))==smap(adjlist(k, 2))
        hasedge(k)=1;
    end
end
p = 0;
p = sum(log(pE(hasedge==1)));
p = p+sum(log(1-pE(hasedge==0)));