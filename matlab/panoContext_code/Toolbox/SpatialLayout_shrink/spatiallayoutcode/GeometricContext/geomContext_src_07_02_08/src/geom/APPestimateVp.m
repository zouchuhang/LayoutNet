function [v2, sigma, p, hpos] = APPestimateVp(lines, imsize, DO_DISPLAY)
% [v2, sigma, p, hpos] = APPestimateVp(lines, imsize, DO_DISPLAY)
% Estimates the principal vanishing points, as in Video Compass [Kosecka
% 2002], except that more (or fewer) than 3 principal vps could be found
% 
% Input:
% lines([x1 x2 y1 y2 angle r])
% imsize: size of image
% DO_DISPLAY (optional): whether to create display figures (default=0)
% Output:
% v2(nvp, [x y]) - the found vanishing points; last is outliers; vp are in
%                  units of pixels/imsize(1), with image upper-left at (0,0)                  
% sigma(nvp) - the variance for each vp
% p(nlines, nvp) - the confidence that each line belongs to each vp
% hpos - horizon position (0 is top of image)
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

if ~exist('DO_DISPLAY')
    DO_DISPLAY = 0;
end

nlines = size(lines, 1);

x1 = [lines(:, [1 3]) ones(size(lines, 1), 1)];
x2 = [lines(:, [2 4]) ones(size(lines, 1), 1)];

% get plane normals for line segments
l = cross(x1, x2);
l = l ./ repmat(sqrt(sum(l.^2,2)), 1, 3);

% make theta range from 0 to pi (instead of -pi to pi)
nbins = 60;
theta = mod(lines(:,5), pi);

% get histogram of thetas
binwidth = pi/nbins;
bincenters = [(binwidth/2):binwidth:(pi-binwidth/2)];
hist_theta = hist(theta, bincenters);

% smooth histogram
for b = 1:nbins
    hist_theta(b) = sum(hist_theta(mod([(b-1):(b+1)]-1, nbins)+1) .* [0.25 0.5 0.25]);
end

% compute curvature of histogram
s = 9;
C = zeros(1, nbins);
for b = 1:nbins
    C(b) = hist_theta(b) - mean(hist_theta(mod([(b-(s-1)/2):(b+(s-1)/2)]-1, nbins)+1));
end


% find dominant peaks surrounded by zero crossings
zc_pos = find((C > 0) & ([C(end) C(1:end-1)]<0));
zc_neg = find((C < 0) & ([C(end) C(1:end-1)]>0));

ngroups = length(zc_pos);
bc = bincenters + pi/nbins/2;
if zc_neg(1) < zc_pos(1)
    i1 = round((zc_pos(end)+zc_neg(end))/2);
    i2 = round((zc_neg(1) + zc_pos(1))/2);
    groups{1} = find((theta > bc(i1)) | (theta < bc(i2)));
    for i = 2:ngroups
        i1 = round((zc_pos(i-1)+zc_neg(i-1))/2);
        i2 = round((zc_neg(i) + zc_pos(i))/2);     
        groups{i} = find((theta > bc(i1)) & (theta < bc(i2)));
    end

else
    for i = 1:ngroups
        i1 = ceil((zc_pos(i)+zc_neg(max(i-1,1)))/2);
        i2 = ceil((zc_neg(i) + zc_pos(min(i+1,ngroups)))/2);
        groups{i} = find((theta > bc(i1)) & (theta < bc(i2)));       
    end
end

% remove groups with few segments 
remove = [];
thresh = max(0.05*nlines, 5);
%thresh = 6;
for i = 1:ngroups
    if length(groups{i}) < thresh
        %disp(num2str([length(groups{i}) nlines]))
        remove(end+1) = i;
    end
end
groups(remove) = [];
ngroups = length(groups);


% initialize EM
sigma = ones(1, ngroups); % sigma = variance here
p = zeros(nlines, ngroups); % p = P(v | l) = P(l|v)P(v)/P(l)
for i = 1:ngroups
    p(groups{i}, i) = 1;
end
%p(:, end) = 1-sum(p(:, 1:end-1), 2);
A = l;
v = [];
sigma = [];
for i = 1:ngroups
    normp = p(:, i) / sum(p(:, i));
    W = diag(normp);
    [eigV, lambda] = eig(A'*W'*W*A);
    [tmp, smallest] = min(diag(lambda));
    v(1:3,i) = eigV(:, smallest);
    sp = sort(normp, 'descend');
    sp = sum(sp(1:min(length(sp), 2)));
    sigma(i) = normp' * (l*v(:,i)).^2 / (1-sum(sp));
end

% create outlier groups
nadd = 3;
ngroups = ngroups + nadd;
sigma(end+(1:nadd)) = 0.2;
tmpv = [0 0 1]';
v(1:3, end+1) = tmpv / sqrt(sum(tmpv.^2));
tmpv = [1 0 1]';
v(1:3, end+1) = tmpv / sqrt(sum(tmpv.^2));
tmpv = [-1 0 1]';
v(1:3, end+1) = tmpv / sqrt(sum(tmpv.^2));
p(:, end+(1:nadd)) = 0;
pv = ones(1, ngroups);

oldv = v;

% do EM
for iter = 1:15
    oldp = p;
    pv = pv / sum(pv);
    % p(l|v) = 1/sqrt(2*pi)/sigma*exp(-(l'v).^2/sigma^2/2) 
    S = repmat(sigma, nlines, 1) + 1E-10;
    plv = exp(-(l*v).^2 ./ S / 2) ./ sqrt(S) + 1E-10;        
    p = plv .* repmat(pv, nlines, 1);    
    p = p ./ repmat(sum(p, 2), 1, ngroups);  
    pv = sum(p, 1);
    nmembers = zeros(ngroups,1);
    for i = 1:ngroups
        normp = p(:, i) / sum(p(:, i));
        W = diag(normp);
        [eigV, lambda] = eig(A'*W'*W*A);
        [tmp, smallest] = min(diag(lambda));
        v(1:3,i) = eigV(:, smallest);            
        sp = sort(normp, 'descend');
        sp = sum(sp(1:min(length(sp), 2)));
        if (1-sum(sp)) > 0
            sigma(i) = normp' * (l*v(:,i)).^2 / (1-sum(sp)); 
        else
            sigma(i) = Inf;                       
        end
    end        

    % remove duplicate groups
    remove =[];
    for i = 1:ngroups        
        for j = i+1:ngroups
            if (v(1:3, i)'*v(1:3,j) > 0.995) || (pv(i)*nlines <= 3) || (sigma(i) > 10) 
                remove(end+1) = i;
                break;
            end
        end
    end
    if length(remove)>0
        p(:, remove) = [];
        sigma(remove) = [];
        v(:, remove) = [];
        pv(remove) = [];
        ngroups = size(p, 2);
    end
    
    % break when v converges
    if all(size(v)==size(oldv)) && min(diag(oldv'*v)) > 0.999
        break;
    end    
    oldv = v;
    
end

% convert vanishing directions to [x y 1] form with (0,0) in upper-left of
% image
%A = inv([1/imsize(1) 0 -1/imsize(2)/2*imsize(1); 0 1/imsize(1) -1/imsize(1)/2*imsize(1); 0 0 1]);
%v2 = (A*v)';
v2 = v';
v2(:, 3) = v2(:, 3)+1E-4;
% sigma = sigma ./ (v2(:, 3).^2)';
v2 = v2 ./ repmat(v2(:, 3), 1, 3);
v2(:, 1) = v2(:, 1) + imsize(2)/imsize(1)/2;
v2(:, 2) = v2(:, 2) + 1/2;

if nargout > 3
    hpos = APPvp2horizon(v2, sigma, p, imsize);
end
%v2 = v;

if DO_DISPLAY

    figure(1), hold off, plot(bincenters, hist_theta, 'b');    
    figure(1), hold on, plot(bincenters, hist_theta, 'r');
    drawnow;

    edge_im = zeros(imsize);
    lines_nnorm(:, [1 2]) = lines(:, [1 2])*imsize(1) + imsize(2)/2;
    lines_nnorm(:, [3 4]) = lines(:, [3 4])*imsize(1) + imsize(1)/2;
    %lines_nnorm = lines;
    for i = 1:length(groups)
        if length(groups{i})>0
            edge_im = draw_line_image2(edge_im, lines_nnorm(groups{i}, 1:4)', i);
        end
    end
    figure(2), imshow(255-label2rgb(edge_im));  
    
    % plot new line memberships
    [tmp, bestv] = max(p, [], 2);    
    edge_im = zeros(imsize);
    for i = 1:ngroups
        groups{i} = find(bestv==i);
        if length(groups{i})>0
            edge_im = draw_line_image2(edge_im, lines_nnorm(groups{i}, 1:4)', i);
        end
        %disp(length(groups{i}))
    end
    figure(3), imshow(255-label2rgb(edge_im));    
    drawnow; 
end


