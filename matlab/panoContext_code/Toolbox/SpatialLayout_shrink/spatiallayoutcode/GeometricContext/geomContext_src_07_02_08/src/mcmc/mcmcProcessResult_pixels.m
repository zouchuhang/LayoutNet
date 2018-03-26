function [vacc, hacc, vcm, hcm] = mcmcProcessResult_pixels(imsegs, pg)
% same as mcmcProcessResult except that pg is a pixel confidence map 

vcm = zeros(3);
hcm = zeros(5);
vtotal = 0;
htotal = 0;
vacc = 0;
hacc = 0;
for f = 1:numel(pg)                        
    
    [imh, imw, nc] = size(pg{f});
    pg{f} = reshape(pg{f}, [imh*imw nc]);
    
    pgv = [pg{f}(:, 1) sum(pg{f}(:, 2:6), 2) pg{f}(:, 7)];
    pgh = pg{f}(:, 2:6);
    
    vlab = imsegs(f).vert_labels(imsegs(f).segimage(:));
    hlab = imsegs(f).horz_labels(imsegs(f).segimage(:));
    npix = numel(vlab);
    vlab = vlab(:);
    hlab = hlab(:);
    
    [maxval, vmax] = max(pgv, [], 2);
    [maxval, hmax] = max(pgh, [], 2);    
    
    vacc = vacc + sum(vlab==vmax)/npix;
    vtotal = vtotal + sum(vlab>0)/npix;
    
    hacc = hacc + sum((hlab==hmax) & (vlab==2))/npix;
    htotal = htotal + sum((hlab>0) & (vlab==2))/npix;        
    
    for s = 1:numel(vmax)
        if vlab(s)~=0
            vcm(vlab(s), vmax(s)) = vcm(vlab(s), vmax(s)) + 1/npix;
        end
        if hlab(s)~=0
            hcm(hlab(s), hmax(s)) = hcm(hlab(s), hmax(s)) + 1/npix;
        end        
    end
end

vacc = vacc / vtotal;
hacc = hacc / max(htotal, 1E-9);

vcm = vcm ./ max(repmat(sum(vcm, 2), 1, size(vcm, 2)), 1E-9);
hcm = hcm ./ max(repmat(sum(hcm, 2), 1, size(hcm, 2)), 1E-9);    
    
%disp(['vacc: ' num2str(vacc) '   hacc: ' num2str(hacc)])
%disp(num2str(sum(ny, 1) / sum(totalv)))