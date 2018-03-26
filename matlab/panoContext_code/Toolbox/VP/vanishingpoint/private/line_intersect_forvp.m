function vp = line_intersect_forvp(line1, line2, p3, p4)
% two syntax: vp = line_intersect_forvp(line1, line2)
%             vp = line_intersect_forvp(p1, p2, p3, p4)
if nargin==2
    p1 = line1.point1; p2 = line1.point2;
    p3 = line2.point1; p4 = line2.point2;
elseif nargin==4
    p1 = line1; p2 = line2;
end

[vp degen] = line_intersect(p1, p2, p3, p4);

if degen==1
%     INFDIST = 1e7; % distance to a far away point...
%     dir = line1.point1 - line1.point2;
%     dir = dir/norm(dir);
%     vp = line1.point1 + dir*INFDIST;
    dir = p1 - p2;
    dir = dir/norm(dir);
    vp = j * dir;
end
