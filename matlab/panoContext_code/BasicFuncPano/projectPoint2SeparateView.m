function [out2D, valid, division] = projectPoint2SeparateView( xyz, vp, fov, sz )
%PROJECTPOINT2SEPARATEVIEW Project a 3D point on panorama to a 2D perspective crops
%   xyz: 3D points on panorama
%   vp: view point of center of the perspective view
%   fov: field of view of the perspective view
%   sz: size of the perspective view, in pixel
%   out2D: corresponding 2D image coordinates
%   valid: check if the projection is valid. Sometimes a 3D point may
%   intersect with a plane in inverse direction, which is invalid.

R = sz/2 / tan(fov/2);
xyz0 = R*vp;
[ ~, out2D, valid, division ] = xyz2view( xyz, xyz0, sz, sz );

end

