function [SAMPLE_ROT_HYP, SAMPLE_COMP] = transformGndRoomsA( objects, para )
%TRANSFORMGNDROOMS Summary of this function goes here
%   Assuming the first one is room

%% add door empty space
% XYZ2POINTLIST = [4 2 6; 4 2 3; 1 2 3; 1 2 6; 
%                  1 5 6; 1 5 3; 4 5 3; 4 5 6];
% I = find([objects.objtype] == 3);
% for i = 1:length(I)
%     xyz = [objects(I(i)).align(1,:) objects(I(i)).align(7,:)];
%     if abs(xyz(1)-xyz(4))<1
%         if xyz(1)<0 % left wall
%             xyz(4) = xyz(1) + 0.2*(xyz(5)-xyz(2));
%         else % right wall
%             xyz(1) = xyz(4) - 0.2*(xyz(5)-xyz(2));
%         end    
%     else
%         if xyz(2)<0 % back wall
%             xyz(5) = xyz(2) + 0.2*(xyz(4)-xyz(1));
%         else % front wall
%             xyz(2) = xyz(5) - 0.2*(xyz(4)-xyz(1));
%         end
%     end
%     objects(I(i)).out_points_w = xyz(XYZ2POINTLIST);
%     x_w = zeros(7,1);
%     x_w(1:3) = xyz(1:3)';
%     x_w(4:6) = xyz(4:6)'-xyz(1:3)';
%     objects(I(i)).x_w = x_w;
%     objects(I(i)) = addBoundingBox(objects(I(i)),'both');
% end

%%
room_obj = objects(1);
room_xyz = [room_obj.align(1,:) room_obj.align(7,:)];
room_x = room_xyz(4) - room_xyz(1);
room_y = room_xyz(5) - room_xyz(2);

floorheight = objects(1).align(1,3) + 0.001;
horzmap = getOccupationMap( objects, floorheight, 1 );
dt_horzmap = bwdist(horzmap);
el_horzmap = dt_horzmap<5;
CC = bwconncomp(~el_horzmap,8);
initCCnumb = CC.NumObjects;
room_area = room_x*room_y;
room_aspect = room_y/room_x;

%% 
[ x_empty_space, y_empty_space, z_empty_space ] = roomEmptyXYZ( objects );
x_empty_space = max(0, x_empty_space-20);
y_empty_space = max(0, y_empty_space-20);

sample_x_min = max(room_x-x_empty_space, para.min_x);
sample_x_max = max(room_x, para.max_x);
sample_y_min = max(room_y-y_empty_space, para.min_y);
sample_y_max = max(room_y, para.max_y);

ROOM_X_SAMPLE = sample_x_max:-10:sample_x_min;
ROOM_Y_SAMPLE = sample_y_max:-10:sample_y_min;

[ grouping_gnd ] = groupingObjectsA( objects, 10, 0, 0);
if grouping_gnd.xfix
    ROOM_X_SAMPLE = room_x;
end
if grouping_gnd.yfix
    ROOM_Y_SAMPLE = room_y;
end

[ROOM_X_VALUE, ROOM_Y_VALUE] = meshgrid(ROOM_X_SAMPLE, ROOM_Y_SAMPLE);
area_trans = (ROOM_X_VALUE.*ROOM_Y_VALUE)/room_area;
aprt_trans = (ROOM_Y_VALUE./ROOM_X_VALUE)/room_aspect;
PRE_VALID = area_trans<2 & area_trans>0.5 & aprt_trans<2 & aprt_trans>0.5;
PRE_VALID = repmat(PRE_VALID, [1 1 5]);

ROOM_VALID = false( size(ROOM_Y_VALUE,1), size(ROOM_X_VALUE,2), 5);
ROOM_HYPOT = cell( size(ROOM_Y_VALUE,1), size(ROOM_X_VALUE,2), 5);

for pushrule = 0:4
    [ grouping_gnd ] = groupingObjectsA( objects, 10, 0, pushrule);
    grouping_int = grouping_gnd;
    previous_x_int = room_x;
    previous_y_int = room_y;
    
%     if grouping_gnd.xfix
%         ROOM_X_SAMPLE = room_x;
%     end
%     if grouping_gnd.yfix
%         ROOM_Y_SAMPLE = room_y;
%     end
    
    for i = 1:length(ROOM_X_SAMPLE)

        sample_x = ROOM_X_SAMPLE(i);
        if sample_x<room_x
            scaleX = sample_x/previous_x_int;
            scaleY = 1;
            objects_int = restoreFromGroupingA( grouping_int, scaleX, scaleY);
            [ grouping_int ] = groupingObjectsA( objects_int, 10, 0, pushrule);
            previous_x_int = sample_x;
            previous_y_int = room_y;
        end

        grouping_cur = grouping_int;
        previous_x_cur = previous_x_int;
        previous_y_cur = previous_y_int;

        for j = 1:length(ROOM_Y_SAMPLE)      
            
            sample_y = ROOM_Y_SAMPLE(j);

            scaleX = sample_x/previous_x_cur;
            scaleY = sample_y/previous_y_cur;
            [ objects_cur ] = restoreFromGroupingA( grouping_cur, scaleX, scaleY);

            ROOM_HYPOT{j,i,pushrule+1} = objects_cur;

            if sample_x<room_x || sample_y<room_y
                valid = validateObjectsRoom(ROOM_HYPOT{j,i,pushrule+1}, grouping_cur, initCCnumb);
                ROOM_VALID(j,i,pushrule+1) = valid;
                if ~valid % if invalid, no need to shrink any more;
                    break; 
                end
            else
                ROOM_VALID(j,i,pushrule+1) = true;
            end

            if sample_y<room_y
                [ grouping_cur ] = groupingObjectsA( objects_cur, 10, 0, pushrule);
                previous_x_cur = sample_x;
                previous_y_cur = sample_y;
            end

        end
    end
end

%% sample and uniform scale change
S = false( length(ROOM_Y_SAMPLE), length(ROOM_X_SAMPLE), 5);
S(1:5:length(ROOM_Y_SAMPLE), 1:5:length(ROOM_X_SAMPLE), :) = true;
S = S & ROOM_VALID & PRE_VALID;
SAMPLE_HYPOT = ROOM_HYPOT(S(:));
uniscale = 2.^([-0.25 0 0.25]);
angle = [0/2*pi 1/2*pi 2/2*pi 3/2*pi];
[gridA, gridS] = meshgrid(angle, uniscale);

SAMPLE_ROT_HYP = cell(length(SAMPLE_HYPOT), numel(gridA));
for sid = 1:length(SAMPLE_HYPOT)
    for tid = 1:numel(gridA)
        ang = gridA(tid);
        scale = gridS(tid);
        hyp = SAMPLE_HYPOT{sid};        
        for oid = 1:length(hyp)
            points = hyp(oid).points;
            hyp(oid).points = points .* scale;

            out_points_w = hyp(oid).out_points_w;
            hyp(oid).out_points_w = out_points_w .* scale;

            x_w = hyp(oid).x_w;
            x_w(1:6) = x_w(1:6) .* scale;
            hyp(oid).x_w = x_w;

            obb = hyp(oid).obb;
            hyp(oid).obb = obb .* scale;

            align = hyp(oid).align;
            hyp(oid).align = align .* scale;
        end
        SAMPLE_ROT_HYP{sid,tid} = rotateObjects(hyp, ang);
    end
end

SAMPLE_INFO = repmat(struct('objtype',[],'feature',[],'objsize',[],'roomctr',[]), numel(SAMPLE_ROT_HYP), 1);
for i = 1:numel(SAMPLE_ROT_HYP)
    SAMPLE_INFO(i).objtype = horzcat(SAMPLE_ROT_HYP{i}.objtype);
    T = reshape(vertcat(SAMPLE_ROT_HYP{i}.align)', 24, []);
    SAMPLE_INFO(i).feature = T;
    SAMPLE_INFO(i).objsize = [sum((T(1:3,:)-T(4:6,:)).^2,1);
                              sum((T(4:6,:)-T(7:9,:)).^2,1);
                              sum((T(13:15,:)-T(1:3,:)).^2,1)];
    SAMPLE_INFO(i).roomctr = (SAMPLE_ROT_HYP{i}(1).align(1,:) + SAMPLE_ROT_HYP{i}(1).align(7,:))/2;
end

SAMPLE_COMP = repmat(struct('objcenter',[],'objsize',[],'objtype',[]), length(SAMPLE_INFO), 1 );
for i = 1:length(SAMPLE_INFO)
    SAMPLE_COMP(i).objcenter = (SAMPLE_INFO(i).feature(1:3,:)+SAMPLE_INFO(i).feature(19:21,:))/2;
    SAMPLE_COMP(i).objcenter = SAMPLE_COMP(i).objcenter - ...
        repmat(SAMPLE_COMP(i).objcenter(:,1), 1, size(SAMPLE_COMP(i).objcenter, 2));
    
    SAMPLE_COMP(i).objsize = SAMPLE_INFO(i).feature(19:21,:)-SAMPLE_INFO(i).feature(1:3,:)+0.1;
    SAMPLE_COMP(i).objtype = SAMPLE_INFO(i).objtype;
    SAMPLE_COMP(i).objaspectXY = SAMPLE_COMP(i).objsize(2,:)./SAMPLE_COMP(i).objsize(1,:);
    SAMPLE_COMP(i).objregionXY = SAMPLE_COMP(i).objsize(2,:).*SAMPLE_COMP(i).objsize(1,:);
    
    SAMPLE_COMP(i).objcenter = single(SAMPLE_COMP(i).objcenter);
    SAMPLE_COMP(i).objsize = single(SAMPLE_COMP(i).objsize);
    SAMPLE_COMP(i).objtype = single(SAMPLE_COMP(i).objtype);
    SAMPLE_COMP(i).objaspectXY = single(SAMPLE_COMP(i).objaspectXY);
    SAMPLE_COMP(i).objregionXY = single(SAMPLE_COMP(i).objregionXY);
end

end

