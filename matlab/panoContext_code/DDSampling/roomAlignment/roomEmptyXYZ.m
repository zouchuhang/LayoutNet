function [ x_empty_space, y_empty_space, z_empty_space, x_span, y_span, z_span ] = roomEmptyXYZ( objects )
%EMPTYXY Summary of this function goes here
%   Assume the first one is room
%%
room_obj = objects(1);
room_ul = room_obj.align(4,1:2);
room_br = room_obj.align(2,1:2);
width = room_br(1) - room_ul(1);
height = room_ul(2)- room_br(2);
scale = 500/width;
swidth = 500;
sheight = round(height*scale);

objects = objects(2:end);
sample_z = zeros(length(objects),2);
for i = 1:length(objects)
    sample_z(i,:) = objects(i).obb(9:10);
end

%% check lower bound
x_empty_space = inf;
y_empty_space = inf;
candi_z = unique([sample_z(:,1)+0.0001; sample_z(:,2)-0.0001]);
for i = 1:length(candi_z)
    cz = candi_z(i);
    inter_obj_id = sample_z(:,1)<cz & sample_z(:,2)>cz;
    
    horzmap = false(sheight, swidth);
    for j = find(inter_obj_id)'
        x = objects(j).obb([1 3 5 7]);
        y = objects(j).obb([2 4 6 8]);
        xs = (x-room_ul(1))*scale;
        xs(abs(xs-swidth)<1) = swidth;
        if std(xs)<1
            if all(xs<1)
                xs(1:2) = 1;
                xs(3:4) = 0;
            else
                xs(1:2) = swidth-1;
                xs(3:4) = swidth;
            end
        end
        
        ys = (room_ul(2)-y)*scale;
        ys(abs(ys-sheight)<1) = sheight;
        if std(ys)<1
            if all(ys<1)
                ys(1:2) = 1;
                ys(3:4) = 0;
            else
                ys(1:2) = sheight-1;
                ys(3:4) = sheight;
            end
        end
        
        mask = poly2mask( xs, ys, sheight, swidth);
        horzmap = horzmap | mask;
    end
    
    xx = min(sum(~horzmap, 2));
    yy = min(sum(~horzmap, 1));
    if xx<x_empty_space
        x_empty_space = xx;
    end
    if yy<y_empty_space
        y_empty_space = yy;
    end  
end

x_empty_space = x_empty_space/scale;
y_empty_space = y_empty_space/scale;

%% check rectangle on wall


%% for z 
z_empty_space = inf;
for i = 1:size(sample_z,1)
    onfloor = abs(sample_z(i,1)-room_obj.obb(9))<10;
    x = objects(i).obb([1 3 5 7]);
    y = objects(i).obb([2 4 6 8]);
    ispolyg = any(max(x)-min(x) < 0.1) | any(max(y)-min(y) < 0.1);
    if onfloor || ~ispolyg
        z_empty_space = min( z_empty_space, room_obj.obb(10)-sample_z(i,2));
    end
end

%%
% x_empty_space = 1 - x_empty_space/width;
% y_empty_space = 1 - y_empty_space/height;
% z_empty_space = 1 - z_empty_space/(room_obj.align(5,3)-room_obj.align(1,3));
x_span = width;
y_span = height;
z_span = (room_obj.align(5,3)-room_obj.align(1,3));
end

