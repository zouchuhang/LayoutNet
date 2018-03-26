function [xo,yo] = intsecl(x1,y1,x2,y2,tol)
% INTSECL Intersection coordinates of two line segments.
%       [XI,YI] = INTSECL(X1,Y1,X2,Y2) where all 4
%	arguments are 2 by N matrices with coordinates
%	of ends of N pairs of line segments (so that 
%       the command such as PLOT(X1,Y1,X2,Y2) will plot 
%       these pairs of lines).
%       Returns 1 by N vectors XI and YI consisting of 
%       coordinates of intersection points of each of N
%	pairs of lines.
%
%	Special cases:
%	When a line segment is degenerate into a point
%	and does not lie on line through the other segment
%	of a pair returns XI=NaN while YI has the following
%	values: 1 - when the first segment in a pair is
%	degenerate, 2 - second segment, 0 - both segments
%	are degenerate.
%	When a pair of line segments is parallel, returns
%	XI = Inf while YI is 1 for coincident lines,
%	0 - for parallel non-coincident ones.
%	INTSECL(X1,Y1,X2,Y2,TOL) also specifies tolerance
%	in detecting coincident points in different line
%	segments.

%  Copyright (c) 1995 by Kirill K. Pankratov
%       kirill@plume.mit.edu
%       04/15/94, 08/14/94, 05/10/95, 08/23/95
 

 % Defaults and parameters .........................
tol_dflt = 0; % Tolerance for coincident points
is_chk = 1;   % Check input arguments

 % Handle input ....................................
if nargin==0, help intsecl, return, end
if nargin<4           % Check if all 4 entered
  error('  Not enough input arguments')
end
if nargin<5, tol = tol_dflt; end
if tol < 0, is_chk = 0; tol = 0; end

 % Check the format of arguments .......
if is_chk
  [x1,y1,x2,y2] = linechk(x1,y1,x2,y2);
end


 % Auxillary
o2 = ones(2,1);
i_pt1 = []; i_pt2 = []; i_pt12 = [];

 % Make first points origins ...........
xo = x1(1,:);
yo = y1(1,:);
x2 = x2-xo(o2,:);
y2 = y2-yo(o2,:);

 % Differences of first segments .......
a = x1(2,:)-x1(1,:);
b = y1(2,:)-y1(1,:);
s = sqrt(a.^2+b.^2);  % Lengths of first segments
i_pt1 = find(~s);
s(i_pt1) = ones(size(i_pt1));
rr = rand(size(i_pt1));
a(i_pt1) = cos(rr);
b(i_pt1) = sin(rr);

 % Normalize by length .................
a = a./s; b = b./s;

 % Rotate coordinates of the second segment
tmp = x2.*a(o2,:)+y2.*b(o2,:);
y2 = -x2.*b(o2,:)+y2.*a(o2,:);
x2 = tmp;

 % Calculate differences in second segments
s = x2(2,:)-x2(1,:);
tmp = y2(2,:)-y2(1,:);
cc = tmp(i_pt1);

 % Find some degenerate cases .......................

 % Find zeros in differences
izy2 = find(~tmp);
tmp(izy2) = ones(size(izy2));

 % Find degenerate and parallel segments
bool = ~s(izy2);
i_par = izy2(~bool);
i_pt2 = izy2(bool);

bool = abs(y2(1,i_pt2))<=tol;
i_pt2_off = i_pt2(~bool);
i_pt2_on = i_pt2(bool);

if length(i_par)~=0
  bool = abs(y2(1,i_par))<=tol;
  i_par_off = i_par(~bool);
  i_par_on = i_par(bool);
end

 % Calculate intercept with rotated x-axis ..........
tmp = s./tmp;   % Slope
tmp = x2(1,:)-y2(1,:).*tmp;


 % Rotate and translate back to original coordinates
xo = tmp.*a+xo;
yo = tmp.*b+yo;

 % Mark special cases ...................................
 % First segments are degenerate to points
if length(i_pt1)~=0
  bool = ~s(i_pt1) & ~cc;
  i_pt12 = i_pt1(bool);
  i_pt1 = i_pt1(~bool);

  bool = abs(tmp(i_pt1))<=tol;
  i_pt1_on = i_pt1(bool);
  i_pt1_off = i_pt1(~bool);

  xo(i_pt1_on) = x1(1,i_pt1_on);
  yo(i_pt1_on) = y1(1,i_pt1_on);

  oo = ones(size(i_pt1_off));
  xo(i_pt1_off) = nan*oo;
  yo(i_pt1_off) = oo;
end

 % Second segments are degenerate to points ...
if length(i_pt2)~=0
  oo = ones(size(i_pt2_off));
  xo(i_pt2_off) = nan*oo;
  yo(i_pt2_off) = 2*oo;
end

 % Both segments are degenerate ...............
if length(i_pt12)~=0
  bool = x1(i_pt12)==xo(i_pt12);
  i_pt12_on = i_pt12(bool);
  i_pt12_off = i_pt12(~bool);

  xo(i_pt12_on) = x1(1,i_pt12_on);
  yo(i_pt12_on) = y1(1,i_pt12_on);

  oo = ones(size(i_pt12_off));
  xo(i_pt12_off) = nan*oo;
  yo(i_pt12_off) = 0*oo;
end

 % Parallel segments .........................
if length(i_par)~=0
  oo = ones(size(i_par_on));
  xo(i_par_on) = inf*oo;
  yo(i_par_on) = oo;

  oo = ones(size(i_par_off));
  xo(i_par_off) = inf*oo;
  yo(i_par_off) = 0*oo;
end




