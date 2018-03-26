function [ valid ] = validateObjectsRoom( objects, grouping, tolccnumb )
%VALIDATEOBJECTSROOM Summary of this function goes here
%   Detailed explanation goes here
% %% add door empty space
% XYZ2POINTLIST = [4 2 6; 4 2 3; 1 2 3; 1 2 6; 
%                  1 5 6; 1 5 3; 4 5 3; 4 5 6];
% I = find([objects.objtype] == 3);
% for i = 1:length(I)
%     xyz = [objects(I(i)).align(1,:) objects(I(i)).align(7,:)];
%     if abs(xyz(1)-xyz(4))<1
%         if xyz(1)<0 % left wall
%             xyz(4) = xyz(1) + 0.6*(xyz(5)-xyz(2));
%         else % right wall
%             xyz(1) = xyz(4) - 0.6*(xyz(5)-xyz(2));
%         end    
%     else
%         if xyz(2)<0 % back wall
%             xyz(5) = xyz(2) + 0.6*(xyz(4)-xyz(1));
%         else % front wall
%             xyz(2) = xyz(5) - 0.6*(xyz(4)-xyz(1));
%         end
%     end
%     objects(I(i)).out_points_w = xyz(XYZ2POINTLIST);
%     x_w = zeros(7,1);
%     x_w(1:3) = xyz(1:3)';
%     x_w(4:6) = xyz(4:6)'-xyz(1:3)';
%     objects(I(i)).x_w = x_w;
%     objects(I(i)) = addBoundingBox(objects(I(i)),'both');
% end

%% check intersection
obbs = obbExpansion(horzcat(objects.obb), 5);
groupid = grouping.groupid;

valid = true;
for i = 2:length(objects)-1
    for j = i+1:length(objects)
        if groupid(i) ~= groupid(j) && ~(objects(i).objtype==3 && objects(j).objtype==3)
            if judgeObbConnect(obbs(:,i), obbs(:,j), 0) 
                valid = false;
                return;
            end
        end
    end
end

%% check out of room
room_min_x = objects(1).align(1,1)-1;
room_max_x = objects(1).align(7,1)+1;
room_min_y = objects(1).align(1,2)-1;
room_max_y = objects(1).align(7,2)+1;

obbs = horzcat(objects.obb)';
bound = [min(obbs(:,[1 3 5 7]),[],2) max(obbs(:,[1 3 5 7]),[],2) min(obbs(:,[2 4 6 8]),[],2) max(obbs(:,[2 4 6 8]),[],2)];
if any(bound(:,1)<room_min_x) || any(bound(:,2)>room_max_x) ...
        || any(bound(:,3)<room_min_y) || any(bound(:,4)>room_max_y)
    valid = false;
    return;
end

%% check connectivity
floorheight = objects(1).align(1,3) + 0.001;
horzmap = getOccupationMap( objects, floorheight, 1 );
dt_horzmap = bwdist(horzmap);
el_horzmap = dt_horzmap<5;

CC = bwconncomp(~el_horzmap,8);
if CC.NumObjects<=tolccnumb
    return;
else
    compsize = zeros(CC.NumObjects,1);
    for i = 1:CC.NumObjects
        compsize(i) = length(CC.PixelIdxList{i});
    end
    if sum(compsize>100)>tolccnumb
        valid = false;
        return;
    end
end

end

