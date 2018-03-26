function r = isinaxes(point, axes)
%ISINAXES determines whether or not a point is in an axes.
%   ISINAXES returns true if POINT is in AXES and false
%   if it is not. POINT is a CurrentPoint. This utility routine
%   works only checks X and Y coordinates.

%   Copyright 2002-2004 The MathWorks, Inc.
%   $Revision: 1.1.4.2 $  $Date: 2004/04/16 22:22:03 $

xr = get(axes,'Xlim');
yr = get(axes,'Ylim');
cx = point(1,1);
cy = point(1,2);
if cx >= xr(1) && cx <= xr(2) && cy >= yr(1) && cy <= yr(2)
    r = true;
else
    r = false;
end;
 
