function [is,in,un] = interval(x1,x2)
% Intersection and union of 2 intervals.
%	[IS,IN,UN] = INTERVAL(X1,X2) calculates pair-wise
%	intersection IN and union UN of N pairs of
%	intervals with coordinates X1 and X2 (both are
%	2 by N vectors). Returns 1 by N boolean vector IS
%	equal to 1 if intervals have non-empty intersection
%	and 0 if they don't.

%  Copyright (c) 1995 by Kirill K. Pankratov,
%       kirill@plume.mit.edu.
%       08/24/95

 % Handle input ...........................
if nargin==0, help interval, return, end
if nargin==1
  un = x1;
else
  un = [x1; x2];
end

[in,un] = sort(un);     % Sort both intervals together
un = un(1:2,:)-1;
is = sum(floor(un/2));  % Check for [0 0 1 1] or [1 1 0 0]
is = (is==1);
ii = find(in(2,:)==in(3,:));
is(ii) = .5*ones(size(ii));

 % Extract intersection and union from sorted coordinates
if nargout>1
  un = in([1 4],:);
  in = in(2:3,:);
  in(:,~is) = flipud(in(:,~is));
end

