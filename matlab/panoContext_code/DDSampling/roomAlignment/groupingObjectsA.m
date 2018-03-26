function [ grouping ] = groupingObjectsA( objects, expand, interthres, pushrule )
%GROUPINGOBJECTS Summary of this function goes here
%   Detailed explanation goes here
if nargin<4
    pushrule = 0;
end

%% add door space
XYZ2POINTLIST = [4 2 6; 4 2 3; 1 2 3; 1 2 6; 
                 1 5 6; 1 5 3; 4 5 3; 4 5 6];
I = find([objects.objtype] == 3);
oriObj = objects(I);
for i = 1:length(I)
    xyz = [objects(I(i)).align(1,:) objects(I(i)).align(7,:)];
    if abs(xyz(1)-xyz(4))<1
        if xyz(1)<0 % left wall
            xyz(4) = xyz(1) + 0.2*(xyz(5)-xyz(2));
        else % right wall
            xyz(1) = xyz(4) - 0.2*(xyz(5)-xyz(2));
        end    
    else
        if xyz(2)<0 % back wall
            xyz(5) = xyz(2) + 0.2*(xyz(4)-xyz(1));
        else % front wall
            xyz(2) = xyz(5) - 0.2*(xyz(4)-xyz(1));
        end
    end
    objects(I(i)).out_points_w = xyz(XYZ2POINTLIST);
    x_w = zeros(7,1);
    x_w(1:3) = xyz(1:3)';
    x_w(4:6) = xyz(4:6)'-xyz(1:3)';
    objects(I(i)).x_w = x_w;
    objects(I(i)) = addBoundingBox(objects(I(i)),'both');
end


%% contacting objects are grouped together
groupid = 1:length(objects);
groupbbox = obbExpansion( horzcat(objects.obb), expand);

for i = 3:length(objects)
    for j = 2:i-1
        if judgeObbConnect(groupbbox(:,i), groupbbox(:,j), interthres)
            groupid(groupid==groupid(j)) = i;
        end
    end
end

XYZ2POINTLIST = [4 2 6; 4 2 3; 1 2 3; 1 2 6; 
                 1 5 6; 1 5 3; 4 5 3; 4 5 6];
groupalign = zeros(length(objects),6);
for i = 1:length(objects)
    groupalign(i,:) = [objects(i).align(1,:) objects(i).align(7,:)];
end
groupobjects = struct([]);
groupbackid = [];
for oid = 1:length(objects)
    if groupid(oid)==oid
        xyz = [min(groupalign(groupid==oid,1:3),[],1) max(groupalign(groupid==oid,4:6),[],1)];        
        obj.out_points_w = xyz(XYZ2POINTLIST);
        obj.type = 3;
        obj.name = 'group';
        obj.objtype = 1;
        groupobjects = [groupobjects; obj];
        groupbackid = [groupbackid; oid];
    end
end
groupobjects = addBoundingBox(groupobjects,'align');
[ rep_point, ~] = objectRepPoint( groupobjects );

objects(I) = oriObj;
offobjects = objects;
for oid = 1:length(objects)
    ref = rep_point(groupbackid==groupid(oid),1:3);
    p = offobjects(oid,:).out_points_w;
    offobjects(oid,:).out_points_w = p - repmat(ref, size(p,1), 1);
    x = offobjects(oid,:).x_w;
    x(1:3) = x(1:3) - ref';
    offobjects(oid,:).x_w = x;
end

grouping.groupobjects = groupobjects;
grouping.groupbackid = groupbackid;
grouping.rep_point = rep_point;
grouping.offobjects = offobjects;
grouping.groupid = groupid;

%% check whether there is a group connecting two facing walls
xfix = false;
yfix = false;
room_xyz = groupobjects(1).align(7,:) - groupobjects(1).align(1,:);
for i = 2:length(groupobjects)
    obj_xyz = groupobjects(i).align(7,:) - groupobjects(i).align(1,:);
    if abs(room_xyz(1)-obj_xyz(1))<1
        xfix = true;
    end
    if abs(room_xyz(2)-obj_xyz(2))<1
        yfix = true;
    end
end
grouping.xfix = xfix;
grouping.yfix = yfix;

%% wall
room_xy = [objects(1).align(1,1:2) objects(1).align(7,1:2)];
J = [3 4 1 2];
K = [1 2 1 2];
if pushrule ~= 0
    i = pushrule;
    adjustvalid = abs(rep_point(:,K(i))-room_xy(i))>1 & abs(rep_point(:,K(i))-room_xy(J(i)))>1;
    adjustvalid(1) = false;
    loc_rep_point = rep_point;
    loc_rep_point(adjustvalid,K(i)) = room_xy(i);
    
    offobjects = objects;
    for oid = 1:length(objects)
        ref = loc_rep_point(groupbackid==groupid(oid),1:3);
        p = offobjects(oid,:).out_points_w;
        offobjects(oid,:).out_points_w = p - repmat(ref, size(p,1), 1);
        x = offobjects(oid,:).x_w;
        x(1:3) = x(1:3) - ref';
        offobjects(oid,:).x_w = x;
    end
    
    grouping.rep_point = loc_rep_point;
    grouping.offobjects = offobjects;    
end



end

