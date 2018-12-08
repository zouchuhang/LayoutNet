function X = seg2poly(s1, P)
% function X = seg2poly(s1, P)
% Check if a line segment s1 intersects with a polygon P.
% INPUTS:
%   s is (2 x 2) where
%     s(:,1) is the first point
%     s(:,2) is the the second point of the segment.
%   P is (2 x n) array, each column is a vertices
% OUTPUT
%   X is (2 x m) array, each column is an intersecting point
%
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original 20-May-2010

% Translate so that first point is origin
a = s1(:,1);
M = bsxfun(@minus, P, a);
b = s1(:,2)-a;
% Check if the points are on the left/right side
x = [b(2) -b(1)]*M;
sx = sign(x);
% x -coordinates has opposite signs
ind = sx(1:end-1).*sx(2:end) <= 0;
if any(ind)
    ind = find(ind);
    % cross point to the y-axis (along the segment)
    x1 = x(ind);
    x2 = x(ind+1);
    d = b.'/(b(1)^2+b(2)^2);
    y1 = d*M(:,ind);
    y2 = d*M(:,ind+1);
    dx = x2-x1;
    % We won't bother with the degenerate case of dx=0 and x1=0
    y = (y1.*x2-y2.*x1)./dx;
    % Check if the cross point is inside the segment
    ind = y>=0 & y<1;
    if any(ind)
        X = bsxfun(@plus, a, b*y(ind));
    else
        X = zeros(2,0);
    end
else
    X = zeros(2,0);
end

end % seg2poly