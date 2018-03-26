function getTwoJunctions(imsegs, adjlist)

[boundmap, perim] = mcmcGetSuperpixelBoundaries(imsegs);
boundmap = boundmap{1};
perim = perim{1};

nsp = imsegs.nseg;
allperim = zeros(nsp, 1);
for s = 1:nsp
    allperim(s) = sum(perim(s, :)) + sum(perim(:, s));
end
    
[ny, nx] = size(imsegs.segimage);

nadj = size(adjlist, 1);

junctions = repmat(struct('type', '', 'sp', [], 'theta', [], 'endpts', []), nadj, 1);

for k = 1:size(nadj, 1)
    s1 = adjlist(k, 1);
    s2 = adjlist(k, 2);
    
    if perim(s1, s2) > 0.8*allperim(s1)
        junctions(k).type = 'enclose';
        junctions(k).sp = [s2 s1]; % s2 encloses s1
    elseif perim(s1, s2) > 0.8*allperim(s2)
        junctions(k).type = 'enclose';
        junctions(k).sp = [s1 s2]; % s1 encloses s2
    elseif perim(s1, s2) > 10    
        ind = boundmap{s1, s2};    
        xvals = mod(ind-1, ny)/(nx-1);
        yvals = floor((ind-1)/ny)/(ny-1);
        
        [p1, p1end, p1error, p2, p2end, p2error] = fit2(xvals, yvals);
        
        if p2error < 0.8*p1error
            junctions(k).type = 'L';
            junctions(k).sp = %CODE FROM HERE
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p1, p1end, p1error, p2, p2end, p2error] = fit2(x, y)

p1 = robustfit(x,y);
tx = polyval(x, p1);
p1error = mean(abs(tx-x));

dy = max(y)-min(y);

dx = max(x)-min(x);

bestp2 = zeros(2,2);
besterror = Inf;
bestt = 0;

if dy > dx % vertical edge            
    [sy, sind] = sort(y);    
    sx = x(sind);   
else
    [sx, sind] = sort(x);    
    sy = y(sind);
end
      
step = max(round(numel(x)/100), 1);
for t = 2:step:numel(x)-1
    p2(1, :) = robustfit(sx(1:t),sy(1:t));
    p2(2, :) = robustfit(sx(t:end),sy(t:end));
    tsy(1:t) = polyval(sx(1:t), p2(1, :));
    tsy(t+1:end) = polyval(sx(t+1:end), p2(2, :));
    terror = mean(abs(tsy-sy));
    if terror < besterror
        bestp2 = p2;
        besterror = terror;
        bestt = t;
    end    
end


p2 = bestp2;
p2error = besterror;        
p2end = [p1end(1, :) ; [sx(bestt) sy(bestt)] ; p1end(2, :)];
   