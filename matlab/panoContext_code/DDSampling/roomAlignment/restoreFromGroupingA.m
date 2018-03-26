function [ objects ] = restoreFromGroupingA(  grouping, scaleX, scaleY  )
%RESTOREFROMGROUPING Summary of this function goes here
%   Detailed explanation goes here

groupbackid = grouping.groupbackid;
groupid = grouping.groupid;

%% data
rep_point = grouping.rep_point;
offobjects = grouping.offobjects;

%% restore big room
scaled_rep_point = [rep_point(:,1)*scaleX rep_point(:,2)*scaleY rep_point(:,3)];
objects = offobjects;

ref = scaled_rep_point(groupbackid==groupid(1),1:3);
% p = offobjects(1).out_points_w + repmat(ref, 8, 1);
% p(:,1) = p(:,1) * scaleX;
% p(:,2) = p(:,2) * scaleY;
p = offobjects(1).out_points_w;
p(:,1) = p(:,1) * scaleX;
p(:,2) = p(:,2) * scaleY;
p = p + repmat(ref, 8, 1);

objects(1).out_points_w = p;

p = offobjects(1).x_w;
p(1) = p(1) * scaleX;
p(2) = p(2) * scaleY;
p(1:3) = p(1:3) + ref';
if abs(p(7)-0)<0.1 || abs(p(7)-pi)<0.1 || abs(p(7)+pi)<0.1
    p(4) = p(4)*scaleX;
    p(5) = p(5)*scaleY;
else
    p(4) = p(4)*scaleY;
    p(5) = p(5)*scaleX;
end
objects(1).x_w = p;

%% restore objects
for oid = 2:length(objects)
    ref = scaled_rep_point(groupbackid==groupid(oid),1:3);
    p = offobjects(oid,:).out_points_w;
    p = p + repmat(ref, size(p,1), 1);
    objects(oid,:).out_points_w = p;
    x = offobjects(oid,:).x_w;
    x(1:3) = x(1:3) + ref';
    objects(oid,:).x_w = x;
end
% all_objects{i} = addBoundingBox(objects,'both');    
objects = addBoundingBox(objects,'both');
end

