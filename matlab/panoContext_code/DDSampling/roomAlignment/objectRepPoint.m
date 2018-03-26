function [ rep_point, objrulelist] = objectRepPoint( objects )
%OBJECTREPPOINT Summary of this function goes here
%   return representing point, and rule of computing point
%   0 means center, 1 means wall, 2 means decide by hypothesis

objrulelist = repmat(struct('usewall',false,'usefloor',false,'free',true), 29, 1);
room_obj = objects(1);
room_xyz = [min(room_obj.align,[],1) max(room_obj.align,[],1)];

rep_point = zeros(length(objects),3);
rep_point(1,:) = (room_xyz(4:6)+room_xyz(1:3))/2;
% objrulelist(29).usecenter = true;

for oid = 2:length(objects)
%     x = objects(oid).obb([1 3 5 7]);
%     y = objects(oid).obb([2 4 6 8]);
    x = objects(oid).align(:,1);
    y = objects(oid).align(:,2);
    objrulelist(objects(oid).objtype).free = false;
    
    ispolyg = any(max(x)-min(x) < 0.1) | any(max(y)-min(y) < 0.1);
    onfloor = abs(objects(oid).align(1,3)-room_obj.align(1,3))<10;
       
    obj_xyz = [min(objects(oid).align,[],1) max(objects(oid).align,[],1)];
    obj_point = (obj_xyz(4:6)+obj_xyz(1:3))/2;
       
    if ~ispolyg || onfloor % use floor
        obj_point(3) = room_obj.align(1,3);
        objrulelist(objects(oid).objtype).usefloor = objrulelist(objects(oid).objtype).usefloor | true;
    end
    
    obj_xy = obj_xyz([1 2 4 5]);
    room_xy = room_xyz([1 2 4 5]);
    I = find( abs( obj_xy-room_xy) < 10);
    if ~isempty(I)
        ID = rem((I-1),2)+1;
%         if sum(ID==1)        
        obj_point(ID) = room_xy(I);
        objrulelist(objects(oid).objtype).usewall = objrulelist(objects(oid).objtype).usewall | true;
    end
    rep_point(oid,:) = obj_point;
end

rep_point = [rep_point vertcat(objects.objtype)];

end

