function [boundmap, perim] = mcmcGetSuperpixelBoundaries(imsegs)
% boundmap{nimages}{nseg, nseg}(npixidx) - boundaries between pairs of segs
% perim{nimages}(nseg, nseg) - number of pixels in boundary
% nimages == 1, boundmap{nseg, nseg}(npixidx), perim(nseg, nseg)

for f = 1:numel(imsegs)
    
    si = imsegs(f).segimage;
    [imh, imw] = size(si);
    
    [gx, gy] = gradient(double(si));
    ind = find((gx.^2 + gy.^2) > 0);
    xi = floor((ind-1) / imh) + 1;
    yi = mod(ind-1, imh)+1;
     
    ti = find((xi < 2) | (xi > imw - 1) | (yi < 2) | (yi > imh - 1));
    ind(ti) = [];
    xi(ti) = [];
    yi(ti) = [];
   
    nseg = imsegs(f).nseg;
    boundmap{f} = cell(nseg, nseg);
    perim{f} = zeros(nseg, nseg);
    npix = zeros(nseg, nseg);
       
    for k = 1:numel(ind)        
        x = xi(k);
        y = yi(k);

        k1 = si(y, x);
        k2 = si(y+1, x);
        if k1 < k2
            npix(k1, k2) = npix(k1, k2) + 1;
        else
            k2 = si(y-1, x);
            if k1 < k2
                npix(k1, k2) = npix(k1, k2) + 1;
            else
                k2 = si(y, x+1);
                if k1 < k2
                    npix(k1, k2) = npix(k1, k2) + 1;
                else
                    k2 = si(y, x-1);
                    if k1 < k2
                        npix(k1, k2) = npix(k1, k2) + 1;
                    end            
                end
            end
        end            
                          
    end
       
    for k1 = 1:nseg
        for k2 = k1+1:nseg
            boundmap{f}{k1, k2} = zeros(npix(k1, k2), 1, 'uint32');
        end
    end    
    npix = zeros(nseg, nseg);    
    
    for k = 1:numel(ind)
        x = xi(k);
        y = yi(k);

        k1 = si(y, x);

        k2 = si(y+1, x);
        if k1 < k2
            npix(k1, k2) = npix(k1, k2) + 1;
            boundmap{f}{k1, k2}(npix(k1, k2)) = ind(k);
        else
            k2 = si(y-1, x);
            if k1 < k2
                npix(k1, k2) = npix(k1, k2) + 1;
                boundmap{f}{k1, k2}(npix(k1, k2)) = ind(k);
            else
                k2 = si(y, x+1);
                if k1 < k2
                    npix(k1, k2) = npix(k1, k2) + 1;
                    boundmap{f}{k1, k2}(npix(k1, k2)) = ind(k);
                else
                    k2 = si(y, x-1);
                    if k1 < k2
                        npix(k1, k2) = npix(k1, k2) + 1;
                        boundmap{f}{k1, k2}(npix(k1, k2)) = ind(k);
                    end            
                end
            end
        end    

    end
           
    perim{f} = npix;

end
    
if f==1
    boundmap = boundmap{1};
    perim = perim{1};
end
