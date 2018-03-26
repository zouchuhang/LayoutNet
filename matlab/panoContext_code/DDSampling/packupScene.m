function [ ALL_SCENE ] = packupScene( ALL_HYPS, CAN_POOL )
%PACKUPSCENE Summary of this function goes here
%   Detailed explanation goes here
global config

[H,W] = size(ALL_HYPS);
ALL_SCENE = cell(H,W);
all_objects = CAN_POOL.sel_hyps.objects;
room_obj.out_points_w = CAN_POOL.room.out_points_w;
room_obj.name = CAN_POOL.room.name;
room_obj.type = CAN_POOL.room.type;
room_obj.points = CAN_POOL.room.points;
room_obj.x_w = CAN_POOL.room.x_w;
room_obj = addBoundingBox(room_obj,'both');
room_obj.objtype = config.roomtypeid; % change this to be room id
for i = 1:H*W
    local = room_obj;
    sel_hyps = ALL_HYPS(i).sel_hyps;
    for j = 1:length(sel_hyps)
        if sel_hyps(j).fixed
            newObj = all_objects(sel_hyps(j).selID);
%             name = get_object_type_livingroom(j); % change this to be function
%             name = get_object_type(j); % change this to be function
            name = feval(config.typefunc, j);
            for k = 1:length(newObj)
                newObj(k).name = name{1};
                newObj(k).objtype = j;
            end
            local = [local; newObj];
        end
    end
    ALL_SCENE{i} = local;
end

end

