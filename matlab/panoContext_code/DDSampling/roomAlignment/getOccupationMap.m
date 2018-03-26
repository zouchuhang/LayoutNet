function [ horzmap ] = getOccupationMap( objects, cz, pixel_per_unit )
%GETOCCUPATIONMAP Summary of this function goes here
%   Detailed explanation goes here
room_obj = objects(1);
room_ul = room_obj.align(4,1:2);
room_br = room_obj.align(2,1:2);
width = room_br(1) - room_ul(1);
height = room_ul(2)- room_br(2);

scale = pixel_per_unit;
sheight = round(height*pixel_per_unit);
swidth = round(width*pixel_per_unit);

objects = objects(2:end);
sample_z = zeros(length(objects),2);
for i = 1:length(objects)
    sample_z(i,:) = objects(i).obb(9:10);
end

% cz = candi_z(i);
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



end

