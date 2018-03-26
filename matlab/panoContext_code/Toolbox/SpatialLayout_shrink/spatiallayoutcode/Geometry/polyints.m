function [xo,yo,ind] = polyints(x1,y1,x2,y2)
% POLYINTS  Intersection of 2 polygons.
%	[XO,YO] = POLYINTS(X1,Y1,X2,Y2) Calculates polygon(s)
%	if intersection of polygons with coordinates X1, Y1 
%	and X2, Y2.
%       The resulting polygon(s) is a set of all points which 
%	belong to both P1 and P2: P = P1 & P2.
%	These polygons must be non-self-intersecting and
%	simply connected.
%
%	If these polygons are not intersecting, returns empty.
%	If intersection consist of several disjoint polygons
%	(for non-convex P1 or P2) output vectors XO, YO consist
%	of concatenated cooddinates of these polygons, 
%	separated by NaN.

%  Copyright (c) 1995 by Kirill K. Pankratov,
%       kirill@plume.mit.edu.
%       06/25/95  
 
if nargin==0, help polyints, return, end
xo=[];yo=[];ind=[];
 % Call POLYBOOL with flag=1
[xo,yo,ind] = polybool(x1,y1,x2,y2,1);
