function [acc, cm, classcount, passthresh] = ...
    msAnalyzeResult(imsegs, labels, pg, DO_SP, ignore, thresh, niter)

nclasses = size(pg{1}, 2);
cm = zeros(nclasses, nclasses);
total = 0;
acc = 0;

passthresh = 0;
pttotal = 0;

if ~exist('thresh', 'var') || isempty(thresh)
    thresh = 0;
end

for f = 1:numel(pg)                        

    if exist('niter', 'var') && ~isempty(niter)     
        pg{f} = mean(pg{f}(:, :, 1:niter), 3);    
    end
    
    lab = labels{f}(:);    
    
    if ~isempty(imsegs)
        npixels = imsegs(f).npixels(:) / sum(imsegs(f).npixels);
    end
    
    pg{f}(:, ignore) = 0;

    ind = false(size(lab));
    for k = 1:numel(ignore)
        ind(lab==ignore(k)) = true;
    end
    lab(ind) = 0;
    [maxval, maxlab] = max(pg{f}, [], 2);   
    
    if DO_SP
    
        pttotal = pttotal + sum(npixels.*(lab>0));
        lab(maxval<thresh) = 0; %ignore labels assigned below threshold confidence
        passthresh = passthresh + sum(npixels.*(lab>0));
        
        acc = acc + sum((lab==maxlab).*npixels);
        total = total + sum(npixels((lab>0)));
        
        for s = 1:numel(maxlab)
            if lab(s)~=0 
                cm(lab(s), maxlab(s)) = cm(lab(s), maxlab(s)) + npixels(s);
            end       
        end
    else        
        
        if ~isempty(imsegs)
            maxlab = maxlab(imsegs(f).segimage);        
            maxval = maxval(imsegs(f).segimage);
        end
        
        pttotal = pttotal + sum((lab>0)/numel(lab));
        lab(maxval<thresh) = 0; %ignore labels assigned below threshold confidence
        passthresh = passthresh + sum((lab>0)/numel(lab));
        
%           maxlab2 = maxlab; maxlab2(maxval<0.5) = 0; 
%         figure(2), hold off, imagesc(maxlab2, [0 21]), axis image
%         figure(3), imagesc(maxlab, [0 21]), axis image        
%         figure(4), imagesc(reshape(lab, size(maxlab)), [0 21]), axis image, colormap jet
%         pause
        
        acc = acc + sum(lab==maxlab(:))/numel(lab);
        total = total + sum(lab(:)>0)/numel(lab);
        for y1 = 1:nclasses
            ind = lab==y1;
            if any(ind)
                for y2 = 1:nclasses
                    ind2 = ind & (maxlab(:)==y2);
                    cm(y1, y2) = cm(y1, y2) + sum(ind2)/numel(lab);
                end
            end
        end
    end
end

acc = acc / total;

classcount = sum(cm, 2)';  
classcount = classcount / sum(classcount);

cm = cm ./ max(repmat(sum(cm, 2), 1, size(cm, 2)), 1E-9);

passthresh = passthresh / pttotal;
    
%disp(['vacc: ' num2str(vacc) '   hacc: ' num2str(hacc)])
%disp(num2str(sum(ny, 1) / sum(totalv)))