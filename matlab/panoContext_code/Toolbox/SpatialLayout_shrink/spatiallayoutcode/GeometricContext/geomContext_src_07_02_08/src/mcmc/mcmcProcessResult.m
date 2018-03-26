function [vacc, hacc, vcm, hcm] = mcmcProcessResult(imsegs, pg, niter)

vcm = zeros(3);
hcm = zeros(5);
vtotal = 0;
htotal = 0;
vacc = 0;
hacc = 0;
for f = 1:numel(pg)                        

    if exist('niter') && ~isempty(niter)     
        pg{f} = mean(pg{f}(:, :, 1:niter), 3);    
    end
    
    pgv = [pg{f}(:, 1) sum(pg{f}(:, 2:6), 2) pg{f}(:, 7)];
    pgh = pg{f}(:, 2:6);
    
    vlab = imsegs(f).vert_labels(:);
    hlab = imsegs(f).horz_labels(:);
    npixels = imsegs(f).npixels(:) / sum(imsegs(f).npixels);
    
    [maxval, vmax] = max(pgv, [], 2);
    [maxval, hmax] = max(pgh, [], 2);    
    
    vacc = vacc + sum((vlab==vmax(:)).*npixels);
    vtotal = vtotal + sum(npixels(find(vlab>0)));
    
    hacc = hacc + sum(((hlab==hmax(:)) & (vlab==2)).*npixels);
    htotal = htotal + sum(npixels((hlab>0) & (vlab==2)));
    
    for s = 1:numel(vmax)
        if vlab(s)~=0
            vcm(vlab(s), vmax(s)) = vcm(vlab(s), vmax(s)) + npixels(s)/numel(imsegs);
        end
        if hlab(s)~=0
            hcm(hlab(s), hmax(s)) = hcm(hlab(s), hmax(s)) + npixels(s)/numel(imsegs);
        end        
    end
end

vacc = vacc / vtotal;
hacc = hacc / max(htotal, 1E-9);

vcm = vcm ./ max(repmat(sum(vcm, 2), 1, size(vcm, 2)), 1E-9);
hcm = hcm ./ max(repmat(sum(hcm, 2), 1, size(hcm, 2)), 1E-9);    
    
%disp(['vacc: ' num2str(vacc) '   hacc: ' num2str(hacc)])
%disp(num2str(sum(ny, 1) / sum(totalv)))