function [ lines ] = lineFromTwoPoint( pt1, pt2 )
%LINEFROMTWOPOINT Generate line segment based on two points on panorama
%   pt1, pt2: two points on panorama
%   lines: 
%       1~3-th dim: normal of the line
%       4-th dim: the projection dimension ID
%       5~6-th dim: the u of line segment endpoints in projection plane
%   use paintParameterLine to visualize

numLine = size( pt1, 1);
lines = zeros(numLine, 6);
n = cross( pt1, pt2, 2);
n = n./repmat( sqrt(sum(n.^2,2)), 1, 3);
lines(:,1:3) = n;

areaXY = abs(sum(n.*repmat([0 0 1], [numLine 1]),2));
areaYZ = abs(sum(n.*repmat([1 0 0], [numLine 1]),2));
areaZX = abs(sum(n.*repmat([0 1 0], [numLine 1]),2));
[~, planeIDs] = max([areaXY areaYZ areaZX], [], 2); % 1:XY 2:YZ 3:ZX
lines(:,4) = planeIDs;

for i = 1:numLine
    uv = xyz2uvN([pt1(i,:); pt2(i,:)], lines(i,4));
    umax = max(uv(:,1))+pi;
    umin = min(uv(:,1))+pi;
    if umax-umin>pi
        lines(i,5:6) = [umax umin]/2/pi;
    else
        lines(i,5:6) = [umin umax]/2/pi;
    end
end

end

