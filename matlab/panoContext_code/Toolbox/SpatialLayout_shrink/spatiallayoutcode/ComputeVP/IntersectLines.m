 %Copyright (c) October,15 2008 by Varsha Hedau, UIUC.  All rights reserved.
function [xs,ys] = IntersectLines(lines1,lines2)

x1s = lines1(:,1); x2s = lines1(:,2);
y1s = lines1(:,3); y2s = lines1(:,4);

x3s = lines2(:,1); x4s = lines2(:,2);
y3s = lines2(:,3); y4s = lines2(:,4);

dets1 = x1s.*y2s - x2s.*y1s;
dets2 = x3s.*y4s - x4s.*y3s;

Nx = dets1.*(x3s-x4s) - dets2.*(x1s-x2s);
Ny = dets1.*(y3s-y4s) - dets2.*(y1s-y2s);
D = (x1s-x2s).*(y3s-y4s) - (x3s-x4s).*(y1s-y2s);

xs = Nx./D;
ys = Ny./D;

return;
