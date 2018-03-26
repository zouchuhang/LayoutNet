function [ out3DNorm, out3DPlane ] = projectPointFromSeparateView(  xy, vp, fov, sz  )
%PROJECTPOINTFROMSEPARATEVIEW Project a point on image to panorama 3D coor.
%   xy: points in 2D image
%   vp: view point of center of the perspective view
%   fov: field of view of the perspective view
%   sz: size of the perspective view
%   out3DNorm: normalized 3D vectors in panorama system
%   out3DPlane: 3D points on cutting planes corresponding to the
%   perspective view.
R = sz/2 / tan(fov/2);
xyz0 = R*vp;

[ out3DNorm, out3DPlane ] = xyzFromView( xy, xyz0, sz, sz );

end

