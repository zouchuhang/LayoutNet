function [lineclass angle] = line_belongto_vp(lines, vp, THRES_THETA)
% find lines that belong to the vanishing point

angle = zeros(1, length(lines));
for k = 1:length(lines)
    angle(k) = anglebetween(lines(k),vp);
end

lineclass = angle < THRES_THETA;

%%
function theta = anglebetween(line, targetpoint)
% get angle between targetpoint & line direction

if ~isreal(targetpoint)
    INFDIST = 1e7; % distance to a far away point...
    targetpoint = INFDIST * targetpoint / i;
end
    
midpoint = (line.point1 + line.point2)/2;

v1 = targetpoint - midpoint;
v2 = line.point2 - midpoint;
theta = 180/pi * real(acos(v1*v2' / norm(v1) / norm(v2)));

if theta>90
	theta = 180 - theta;
end
