function [is,S] = iscross(x1,y1,x2,y2,tol)
% ISCROSS  Finds whether pairs of lines cross each other
%	[IS,S] = ISCROSS(X1,Y1,X2,Y2) where arguments X1, Y1,
%	X2, Y2 are all 2 by N matrices are coordinates of
%	ends of the pairs of line segments.
%	Returns vector  IS (1 by N)  consisting of ones if 
%	corresponding pairs cross each other, zeros if they 
%	don't and .5 if an end of one line segment lies on
%	another segment.
%	Also returns a matrix  S (4 by N) with each row 
%	consisting of cross products (double areas of 
%	corresponding triangles) built on the following points:
%	(X2(1,:),Y2(1,:)),(X1(1,:),Y1(1,:)),(X2(2,:),Y2(2,:)),
%	(X2(1,:),Y2(1,:)),(X1(2,:),Y1(2,:)),(X2(2,:),Y2(2,:))
%	(X1(1,:),Y1(1,:)),(X2(1,:),Y2(1,:)),(X1(2,:),Y1(2,:))
%	(X1(1,:),Y1(1,:)),(X2(2,:),Y2(2,:)),(X1(2,:),Y1(2,:))
%	The signs of these 4 areas can be used to determine
%	whether these lines and their continuations cross each
%	other.
%	[IS,S] = ISCROSS(X1,Y1,X2,Y2,TOL) uses tolerance TOL
%	for detecting the crossings (default is 0).

%  Copyright (c) 1995 by Kirill K. Pankratov
%       kirill@plume.mit.edu
%       08/14/94, 05/18/95, 08/25/95

 % Defaults and parameters .......................
tol_dflt = 0; % Tolerance for area calculation
is_chk = 1;   % Check input arguments

 % Handle input ..................................
if nargin==0, help iscross, return, end
if nargin<4           % Check if all 4 entered
  error('  Not enough input arguments')
end
if nargin<5, tol = tol_dflt; end
if tol < 0, is_chk = 0; tol = 0; end

 % Check the format of arguments .................
if is_chk
  [x1,y1,x2,y2] = linechk(x1,y1,x2,y2);
end

len = size(x1,2);
o2 = ones(2,1);

 % Find if the ranges of pairs of segments intersect
[isx,S,A] = interval(x1,x2);
scx = diff(A);
[isy,S,A] = interval(y1,y2);
scy = diff(A);
is = isx & isy;

 % If S values are not needed, extract only those pairs
 % which have intersecting ranges ..............
if nargout < 2
  isx = find(is);  % Indices of pairs to be checked
                   % further
  x1 = x1(:,isx);
  x2 = x2(:,isx);
  y1 = y1(:,isx);
  y2 = y2(:,isx);
  is = is(isx);
  if length(is)==0, is = zeros(1,len); return, end
  scx = scx(isx);
  scy = scy(isx);
end

 % Rescale by ranges ...........................
x1 = x1.*scx(o2,:);
x2 = x2.*scx(o2,:);
y1 = y1.*scy(o2,:);
y2 = y2.*scy(o2,:);


 % Calculate areas .............................
S = zeros(4,length(scx));
S(1,:) = (x2(1,:)-x1(1,:)).*(y2(2,:)-y1(1,:));
S(1,:) = S(1,:)-(x2(2,:)-x1(1,:)).*(y2(1,:)-y1(1,:));

S(2,:) = (x2(1,:)-x1(2,:)).*(y2(2,:)-y1(2,:));
S(2,:) = S(2,:)-(x2(2,:)-x1(2,:)).*(y2(1,:)-y1(2,:));

S(3,:) = (x1(1,:)-x2(1,:)).*(y1(2,:)-y2(1,:));
S(3,:) = S(3,:)-(x1(2,:)-x2(1,:)).*(y1(1,:)-y2(1,:));

S(4,:) = (x1(1,:)-x2(2,:)).*(y1(2,:)-y2(2,:));
S(4,:) = S(4,:)-(x1(2,:)-x2(2,:)).*(y1(1,:)-y2(2,:));


 % Find if they cross each other ...............
is = (S(1,:).*S(2,:)<=0)&(S(3,:).*S(4,:)<=0);


 % Find very close to intersection
isy = min(abs(S));
ii = find(isy<=tol & is);
is(ii) = .5*ones(size(ii));

 % Output
if nargout < 2
  isy = zeros(1,len);
  isy(isx) = is;
  is = isy;

else
  isy = scx.*scy;
  ii = find(~isy);
  isy(ii) = ones(size(ii));
  S = S./isy(ones(4,1),:);

end

