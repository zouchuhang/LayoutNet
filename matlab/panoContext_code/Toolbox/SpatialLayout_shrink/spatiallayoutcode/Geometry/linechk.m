function [x1,y1,x2,y2] = linechk(x1,y1,x2,y2)
% LINECHK Input checking for line segments.

%  Copyright (c) 1995 by Kirill K. Pankratov
%       kirill@plume.mit.edu
%       08/22/95,

 % String for transposing
str = ['x1=x1'';'; 'y1=y1'';'; 'x2=x2'';'; 'y2=y2'';'];

 % Sizes
sz = [size(x1); size(y1); size(x2); size(y2)]';
psz = prod(sz);

 % Check x1, y1
if psz(1)~=psz(2)
  error('  Arguments  x1 and y1 must have the same size')
end

 % Check x2, y2
if psz(3)~=psz(3)
  error('  Arguments  x2 and y2 must have the same size')
end

 % Check if any arguments are less than 2 by 1
if any(max(sz)<2)
  error('  Arguments  x1, y1, x2, y2 must be at least 2 by 1 vectors')
end

 % Check if no size is equal to 2
if any(all(sz~=2))
  error('  Arguments  x1, y1, x2, y2 must be 2 by 1 vectors')
end

 % Find aruments to be transposed .............................
ii = find(sz(1,:)~=2);
for jj = 1:length(ii)
  eval(str(ii(jj),:)); % Transpose if neccessary
end
sz(:,ii) = flipud(sz(:,ii));

 % If vectors, extend to 2 by n matrices ......................
n = max(sz(2,:));
on = ones(1,n);
if sz(2,1)<n, x1 = x1(:,on); end
if sz(2,2)<n, y1 = y1(:,on); end
if sz(2,3)<n, x2 = x2(:,on); end
if sz(2,4)<n, y2 = y2(:,on); end

