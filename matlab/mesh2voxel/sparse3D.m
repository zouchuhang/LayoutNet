function S = sparse3D(M)
% SPARSE3D  Converts a 3D binary matrix to a Nx3 coodrinate matrix.
%    S = SPARSE3D(M) Converts the 3D binary matrix (occupancy matrix) M to
%       a Nx3 matrix. M is a 3D matrix where 1-valued cells are considered
%       as full and 0-filled cells are empty. The returned matrix S
%       contains the XYZ coordinates of the cells which are occupied.
% 
%    Example 1
%    ---------
%    The following example creates the occupancy matrix of a very small
%    sphere and then converts it to sparse representation:
% 
%       a = zeros(5,5,5);
%       for i=1:5; x = i-3;
%         for j=1:5; y = j-3;
%           for k=1:5; z = k-3;
%             if sum([ x y z ].^2)<=4; a(i,j,k) = 1; end;
%           end
%         end
%       end
%       sp = sparse3D(a);
%       disp(a);
%       disp(sp);
%
  S = [];
  for z=1:size(M,3)
    [ x y ] = find(M(:,:,z)==1);
    S = vertcat(S,[ x y ones(size(x,1),1)*z ]);
  end
