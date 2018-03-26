function [ score ] = findNearbyObject( all_xyz, obj_xyz )
%FINDNEARBYOBJECT Summary of this function goes here
%   Detailed explanation goes here
N = size(all_xyz,1);

lowxyz = max(all_xyz(:,1:3), repmat(obj_xyz(1:3),N,1));
topxyz = min(all_xyz(:,4:6), repmat(obj_xyz(4:6),N,1));
inter = prod(max(topxyz-lowxyz, 0.01), 2);
area1 = prod(max(all_xyz(:,4:6)-all_xyz(:,1:3), 0.01), 2);
area2 = repmat(prod(max(obj_xyz(4:6)-obj_xyz(1:3), 0.01), 2), N, 1);

score = inter./(area1+area2-inter);
valid = any(topxyz-lowxyz<=0,2);
score(valid) = 0;

end

