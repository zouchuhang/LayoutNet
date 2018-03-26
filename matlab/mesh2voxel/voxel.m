function varargout = voxel(i,d,c,alpha)
% VOXEL Draw a 3-D voxel.
%    VOXEL(S) Draws a voxel centered at the specified point S. S is a
%       coordinate point in the form [ x y z ].
%    VOXEL(S,EDGE) Draws a voxel using the specified EDGE size. If no
%       EDGE is provided, it is set by default to 1. EDGE is a three
%       element vector [dx,dy,dz].
%    VOXEL(S,EDGE,C) Uses the specified colour C to draw the faces of the
%       voxel. C is a character string to specify color (type 'help plot'
%       to see list of valid colors). If no C is provided, it is set by
%       default to blue.
%    VOXEL(S,EDGE,C,ALPHA) Uses the given ALPHA to define the transparency
%       level (1 for opaque, 0 for transparent). If no ALPHA is given, it
%       is set to 1.
%    H = VOXEL(...) Return a vector of handles for the voxel drawn.
%
%
% See also PLOT, PATCH

%   Suresh Joel Apr 15,2003
%           Updated Feb 25, 2004

switch(nargin),
case 0
    disp('Too few arguements for voxel');
    return;
case 1
    d=[ 1 1 1 ]; %default length of side of voxel is 1
    c='b';       %default color of voxel is blue
    alpha=1;
case 2
    c='b';
    alpha=1;
case 3
    alpha=1;
end;

x=[i(1)+[0 0 0 0 d(1) d(1) d(1) d(1)]; ...
        i(2)+[0 0 d(2) d(2) 0 0 d(2) d(2)]; ...
        i(3)+[0 d(3) 0 d(3) 0 d(3) 0 d(3)]]';
h = [];
for n=1:3,
    if n==3,
        x=sortrows(x,[n,1]);
    else
        x=sortrows(x,[n n+1]);
    end;
    temp=x(3,:);
    x(3,:)=x(4,:);
    x(4,:)=temp;
    h1=patch(x(1:4,1),x(1:4,2),x(1:4,3),c);
    set(h1,'FaceAlpha',alpha);
    h = vertcat(h1,h);
    temp=x(7,:);
    x(7,:)=x(8,:);
    x(8,:)=temp;
    h1=patch(x(5:8,1),x(5:8,2),x(5:8,3),c);
    set(h1,'FaceAlpha',alpha);
    h = vertcat(h1,h);
end;

if nargout>0
  varargout{1} = h;
end
