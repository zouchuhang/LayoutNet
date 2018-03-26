function [pt degen] = line_intersect(p1, p2, p3, p4)
% function pt = line_intersect(p1, p2, p3, p4)
% intersection point of line defined by p1 & p2 and line defined by p3 & p4
% http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
x1 = p1(1); y1 = p1(2);
x2 = p2(1); y2 = p2(2);
x3 = p3(1); y3 = p3(2);
x4 = p4(1); y4 = p4(2);

if (y4-y3)*(x2-x1)-(x4-x3)*(y2-y1) == 0
% 	warning('line_intersect.m degenerate --dclee');
% 	pt = (p1 + p2 + p3 + p4)/4;
    pt = [];
    degen = 1;
	return;
end

pt = [ ...
	x1 + (x2-x1) * ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/((y4-y3)*(x2-x1)-(x4-x3)*(y2-y1)) ...
	y1 + (y2-y1) * ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/((y4-y3)*(x2-x1)-(x4-x3)*(y2-y1)) ];
degen = 0;
