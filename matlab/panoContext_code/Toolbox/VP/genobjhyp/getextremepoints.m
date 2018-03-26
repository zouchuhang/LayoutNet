function [pleft pright pbot2 ptop2 pbot3 ptop3 pleft3 pright3] = ...
    getextremepoints(region, vp)
% region: [n x 2] [x y] position of n points

%% get extreme points in image
[point1 point2] = tangentpoint(region, vp{1});
if point1(1) > point2(1)
    pright = point1; pleft = point2;
else
    pright = point2; pleft = point1;
end

[point1 point2] = tangentpoint(region, vp{2});
if point1(2) > point2(2)
    pbot2 = point1; ptop2 = point2;
else
    pbot2 = point2; ptop2 = point1;
end

[point1 point2] = tangentpoint(region, vp{3});
if point1(2) > point2(2)
    pbot3 = point1; ptop3 = point2;
else
    pbot3 = point2; ptop3 = point1;
end
if point1(1) > point2(1)
    pright3 = point1; pleft3 = point2;
else
    pright3 = point2; pleft3 = point1;
end

function [point1 point2] = tangentpoint(region, pivot)
% find tangent point of a region
% region: [n x 2] [x y] position of n points
% pivot: [1 x 2] [x y] position of pivot

n = size(region,1);

xslope = region(:,1) - repmat(pivot(1), n, 1);
yslope = region(:,2) - repmat(pivot(2), n, 1);

ang = atan2(yslope, xslope);

% mean angle in circular statistics
% (http://en.wikipedia.org/wiki/Mean_of_circular_quantities)
meanang = atan2(mean(sin(ang)), mean(cos(ang)));

% circular distance
cd = circ_dist(ang, meanang);

[c maxidx] = max(cd);
[c minidx] = min(cd);

point1 = region(maxidx, :);
point2 = region(minidx, :);


%%
function r =  circ_dist(x,y)
%
% r = circ_dist(alpha, beta)
%   Pairwise difference x_i-y_i around the circle computed efficiently.
%
%   Input:
%     alpha      sample of linear random variable
%     beta       sample of linear random variable or one single angle
%
%   Output:
%     r       matrix with differences
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 3/19/2009
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


if size(x,1)~=size(y,1) && size(x,2)~=size(y,2) && length(y)~=1
  error('Input dimensions do not match.')
end

r = angle(exp(1i*x)./exp(1i*y));
