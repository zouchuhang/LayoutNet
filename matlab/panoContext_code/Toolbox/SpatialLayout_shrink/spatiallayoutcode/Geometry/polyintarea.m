function intarea = polyintarea(x1,y1,x2,y2,flag)
% POLYINTAREA  Intersection area of 2 polygons.
%	INTAREA = POLYINTS(X1,Y1,X2,Y2) Calculates polygon(s)
%	if intersection of polygons with coordinates X1, Y1 
%	and X2, Y2.
%   Code computes intersection polygon similar to polyints and
%   and uses it to compute intersection area

%  Copyright (c) 1995 by Kirill K. Pankratov,
%       kirill@plume.mit.edu.
%       06/25/95  
 
if nargin==0, help polyints, return, end
xo=[];yo=[];ind=[];
 % Call POLYBOOL with flag=1
[xo,yo,ind] = polybool(x1,y1,x2,y2,1,flag);
if numel(xo)==0
    intarea=0;
    return;
end

ii = find(isnan(xo));
if numel(ii)==0
    intarea=polyarea([xo;xo(1)],[yo;yo(1)]);
    return;
end

intarea=0;
ii=[0;ii;length(xo)+1];
for i=1:numel(ii)-1
    seg_xo=xo(ii(i)+1:ii(i+1)-1);
    seg_yo=yo(ii(i)+1:ii(i+1)-1);
    if numel(seg_xo)>0
        intarea = intarea+polyarea([seg_xo;seg_xo(1)],[seg_yo;seg_yo(1)]);
    end
end

return;
