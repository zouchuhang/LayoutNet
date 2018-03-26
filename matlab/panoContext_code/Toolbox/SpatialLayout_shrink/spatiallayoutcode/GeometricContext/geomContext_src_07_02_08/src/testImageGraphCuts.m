function [labv, labh] = testImageGraphCuts(spdata, edata, adjlist, vclassifierSP, hclassifierSP, eclassifier, ecal)

% probability of superpixel main labels
pvSP = test_boosted_dt_mc(vclassifierSP, spdata);
pvSP = 1 ./ (1+exp(-pvSP));
pvSP = pvSP ./ repmat(sum(pvSP, 2), 1, size(pvSP, 2));
[tmp, vmax] = max(pvSP, [], 2);

% probability of superpixel sub labels
phSP = test_boosted_dt_mc(hclassifierSP, spdata);
phSP = 1 ./ (1+exp(-phSP));
phSP = phSP ./ repmat(sum(phSP, 2), 1, size(phSP, 2));
[tmp, hmax] = max(phSP, [], 2);

% edge probability
pE = test_boosted_dt_mc(eclassifier, edata);
pE = 1 ./ (1+exp(ecal(1)*pE+ecal(2)));

%labv = alphaExpansion(pvSP, pE, adjlist);
labh = hmax;
labv = vmax;
% vind = find(labv==2);
% vaind = find((labv(adjlist(:, 1)) == 2) & (labv(adjlist(:, 2)) == 2));
% adjlisth = adjlist(vaind, :);
% for k = 1:size(adjlisth, 1)
%     s1 = find(adjlisth(k, 1)==vind);
%     s2 = find(adjlisth(k, 2)==vind);
%     adjlisth(k, :) = [s1 s2];
% end
% pEh = pE(vaind);
% phSP = phSP(vind, :);
% 
% labh = hmax;
% tmplabh = alphaExpansion(phSP, pEh, adjlisth);   
% labh(vind) = tmplabh;

%labh = alphaExpansion(phSP, pE, adjlist);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function lab = alphaExpansion(plab, pE, adjlist)

nsp = size(plab, 1);
nadj = size(adjlist, 1);

epenalty = max(log(pE)-log(1-pE), 0);% -log(1-pE); % 

[tmp, lab] = max(plab, [], 2);
success = 1;
while success
    success = 0;
    for y = 1:3
        aux = find(lab(adjlist(:,1))~=lab(adjlist(:, 2)));
        naux = numel(aux);
        nnodes = nsp + naux;
        
        STE = zeros(2, nnodes);
        STE(1, 1:nsp) = -log(plab(:, y)');
        aind = find(lab==y);
        naind = find(lab~=y);
        STE(2, aind) = Inf;
        STE(2, naind) = -log(plab(nsp*(lab(naind)-1)+naind)); % -log P(current label)
        STE(2, nsp+1:end) = epenalty(aux);

        % nodes include main nodes and aux nodes for dealing with edges
        % between different-label nodes
        E_n = zeros(2, nadj+naux);
        E_n(1:2, 1:nadj) = adjlist';        
        E_n(2, (nadj+1):(nadj+naux)) = E_n(2, aux);
        E_n(2, aux) = nsp+[1:naux];
        E_n(1, (nadj+1):(nadj+naux)) = E_n(2, aux);
        
        E_w = zeros(2, nadj+naux);
        for k = setdiff([1:nadj], aux')            
            if lab(adjlist(k, 1)==y)
                E_w(:, k) = 0;
            else
                E_w(:, k) = epenalty(k);
            end
        end        
        for k = 1:numel(aux)
            k1 = aux(k);
            k2 = k+nadj;
            if lab(adjlist(k1, 1)~=y)
                E_w(:, k) = epenalty(k1);
            end
            if lab(adjlist(k1, 2)~=y)
                E_w(:, k2) = epenalty(k1);
            end            
        end
       
        [cut, flow] = vgg_graph_maxflow(uint32(nnodes), int16(100*STE), uint32(E_n), int16(100*E_w));
        
        lab2 = lab;
        lab2(logical(cut(1:nsp))) = y;
        
        energy1 = sum(-log(plab(nsp*(lab-1)+[1:nsp]'))) + sum(epenalty(aux));
        
        aux2 = find(lab2(adjlist(:,1))~=lab2(adjlist(:, 2)));
        energy2 = sum(-log(plab(nsp*(lab2-1)+[1:nsp]'))) + sum(epenalty(aux2));
                        
        if energy2 < energy1
            %disp(['Energy 1 = ' num2str(energy1) '   Energy 2 = ' num2str(energy2)]);
            lab = lab2;
            success = 1;
        end
    end
end
        
 