function [pos, tilt, conf, y, py] = APPvp2horizon(v, vars, p, imsize)
% [pos, tilt, conf, y, py] = APPvp2horizon(v, vars, p, imsize)
%
% Estimates the horizon position based on the vanishing points. Doesn't
% work well.
%
% Input: 
% v - vanishing points from estimateVP
% vars - variances on vanishing point estimates
% p - memberships of lines in vp
% imsize - the size of the grayscale image
% Output:
% pos - y position (0 is top of image) of horizon, NaN if cannot be found
% tilt - the angle of the horizon line, NaN if not set
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

% find vanishing points that occur within image height (or close)
ind = find([v(:, 2) > -0.25 & v(:, 2) < 1.25]);

% do not consider vp with variance greater than varthresh
varthresh = 0.005;
ind(find(vars(ind) > varthresh)) = [];

% do not consider vp with fewer than 6 lines (expected)
memberthresh = 6;
nmembers = sum(p(:, ind), 1);
ind(find(nmembers < memberthresh)) = [];

conf = length(ind);

% if more than two vp are left, choose two with least variance
%[tmp, ordering] = sort(vars(ind), 'ascend');
%ind = ind(ordering(1:min(2, length(ordering))));

v = v(ind, :);
vars = vars(ind);

% compute horizon position (if possible) and tilt (if possible)
if length(ind)==0
    pos = NaN;
    tilt = NaN;
    conf = NaN;
elseif length(ind)>0

    % change variables so variance is along y axis only
    v(:, 1) = v(:, 1) - 0.5*imsize(2)/imsize(1);
    v(:, 2) = v(:, 2) - 0.5;
    for i = 1:length(ind)
        vars(i) = vars(i)*v(i,2)^2 / (v(i,1)^2+v(i,2)^2);    
    end    
    pos = sum(v(:, 2) ./ vars') / sum(1./vars');
%    pos = mean(v(:, 2));
    
%     pos = mean(v(:, 2));
%     conf = std(v(:, 2)) * sqrt(length(ind) / (length(ind)-1));
%     tilt = 0;

%     disp(['vp: ' num2str(length(ind))])
%     if length(ind)>2
%         [p, stats] = robustfit(v(:,1),v(:,2));
%         tilt = p(1);
%     else
%         [p, stats] = robustfit(v(:, 1), v(:, 2));         
%         tilt = p(1);
%     end
%     conf = stats.se(2);
%     pos = p(2);
    %[pos, conf] = polyval(p, 0, S);
    %pos = p(2);
%    conf = mu(2);
    %tilt = p(1);
%     disp(num2str([pos+0.5 tilt conf]));
        
    % get distribution of horizon surrounding most likely position
    if nargout > 3
        y = [pos-0.25:0.005:pos+0.25];
        py = zeros(size(y));
        for i = 1:size(v, 1)
            py = py + log(1 / sqrt(2*pi*vars(i))*exp(-1/2/vars(i)*(y-v(i,2)).^2));
        end
        py = exp(py);
        py = py / sum(py);

        y = y + 0.5;

        % remove sampling points with very low probability
        ind = find(py> (1E-4)/length(py));
        y = y(ind);
        py = py(ind);
        py = py / sum(py);
    end
    
    pos = pos + 0.5;
    %tilt = NaN;    
    
    %figure(1), plot(y, py, 'r')    
end


% elseif length(ind)==1
%     pos = 1-v(1, 2);
%     tilt = NaN;
% end
% else % length(ind)==2
%    midx = (imsize(2)/2/imsize(1));
%    m = (v(1,2)-v(2,2)) / (v(1,1)-v(2,1));
%    b = v(1,2) - m*v(1,1);
%    pos = midx*m+b;
%    tilt = atan(-m);
% end
