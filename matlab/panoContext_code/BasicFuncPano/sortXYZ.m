function [ s_xyz, I ] = sortXYZ( xyz )
%SORTXYZ Sort 3D point in clockwise order
%   The definition of clockwise order for 3D points is ill-posed. To make it well defined,
%   we must first define a viewing direction. This code treats the center
%   of the set of points as the viewing direction (from original point).
ori_xyz = xyz;
xyz = xyz./repmat(sum(xyz.^2,2),1,3);
center = mean(xyz,1);
center = center./norm(center);
% set up z axis at center
z = center;
x = [-center(2) center(1) 0]; x = x./norm(x);
y = cross(z,x,2);
R = diag([1 1 1])/[x' y' z'];

newXYZ = R*(xyz');
A = atan2(newXYZ(2,:), newXYZ(1,:));
[~,I] = sort(A,'ascend');

s_xyz = ori_xyz(I,:);


end

