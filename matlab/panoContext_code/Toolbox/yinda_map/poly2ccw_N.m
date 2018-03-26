function [x2, y2] = poly2ccw_N(x1, y1)
%POLY2CCW_N Convert polygon contour to counterclockwise vertex ordering
%
%   [X2, Y2] = POLY2CCW(X1, Y1) arranges the vertices in the polygonal
%   contour (X1, Y1) in counterclockwise order, returning the result in X2
%   and Y2. If X1 and Y1 can contain multiple contours, represented either
%   as NaN-separated vectors or as cell arrays, then each contour is
%   converted to clockwise ordering.  X2 and Y2 have the same format
%   (NaN-separated vectors or cell arrays) as X1 and Y1.
%
%   Example
%   -------
%   Convert a clockwise-ordered square to counterclockwise ordering.
%
%       x1 = [0 0 1 1 0];
%       y1 = [0 1 1 0 0];
%       ispolycw(x1, y1)
%       [x2, y2] = poly2ccw(x1, y1);
%       ispolycw(x2, y2)
%
%   See also ISPOLYCW, POLY2CW, POLYBOOL.

% Copyright 2004-2009 The MathWorks, Inc.
% $Revision: 1.1.4.3 $  $Date: 2009/08/11 15:44:33 $


input_is_cell = iscell(x1);
if ~input_is_cell
   checkxy_N(x1, y1, mfilename, 'X1', 'Y1', 1, 2)
   input_is_row = (size(x1, 1) == 1);
   [x1, y1] = polysplit_N(x1, y1);
end

x2 = x1;
y2 = y1;

for k = 1:numel(x1)
   if ispolycw_N(x2{k}, y2{k})
      x2{k} = x2{k}(end:-1:1);
      y2{k} = y2{k}(end:-1:1);
   end
end

if ~input_is_cell
   [x2, y2] = polyjoin_N(x2, y2);
   if input_is_row
      x2 = x2';
      y2 = y2';
   end
end
