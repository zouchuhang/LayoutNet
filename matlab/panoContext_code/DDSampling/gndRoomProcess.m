function [ ALL_GNDS_ROT, ALL_GNDS_EMPTYSPACE, ALL_GNDS_REPPOINT, ALL_GNDS_RULES] = gndRoomProcess( folderName )
%GNDROOMPROCESS Summary of this function goes here
%   Detailed explanation goes here
load([folderName '/IMGLIST.mat']);
load([folderName '/ANNO_ALL.mat']);

%% get bounding box
anno_valid = false( 1, length(ANNO_ALL));
for i = 1:length(anno_valid)
    anno_valid(i) = ~isempty(ANNO_ALL(i).ANNO3D) && ANNO_ALL(i).ANNO3D.b_singleroom;
end
valid_anno_id = find(anno_valid);
ALL_GNDS = cell(length(valid_anno_id),1);
for i = 1:length(valid_anno_id)
    id = valid_anno_id(i);
    fprintf('%d\n',id);
    contain3Dcoords = ANNO_ALL(id).ANNO3D.contain3Dcoords;
    gnd_annotation = ANNO_ALL(id).ANNO3D.objects3D(contain3Dcoords);
    x_w = gnd_annotation(1).x_w;
    for j = 1:length(gnd_annotation)
        gnd_annotation(j).objtype = get_object_type_livingroom({gnd_annotation(j).name});
        gnd_annotation(j).x_w = x_w(j*7-6:j*7);
    end
    ALL_GNDS{id} = rmfield( gnd_annotation, ...
        {'out_points_o','x_o','out_points_1','x_1','ori_polygon', ...
        'out_points_2','x_2','points_all','points_adj','out_points_3','x_3'});
%     ALL_GNDS{id} = addBoundingBox(ALL_GNDS{id},'both');
end

%% rotate ground truth
angle = [0 pi/2 pi 3/2*pi];
ALL_GNDS_ROT = cell(length(ALL_GNDS), 4);
for aid = 1:length(ALL_GNDS)
    if isempty(ALL_GNDS{aid})
        continue;
    end
    for rid = 1:1%4
        ALL_GNDS_ROT{aid,rid} = rotateObjects(ALL_GNDS{aid}, angle(rid));
    end
end

%% compute xy scales
ALL_GNDS_EMPTYSPACE = cell(length(ALL_GNDS), 4);

for aid = 1:length(ALL_GNDS)
    if isempty(ALL_GNDS{aid})
        continue;
    end
    for rid = 1:1%4
        [x y z xl yl zl] = roomEmptyXYZ(ALL_GNDS_ROT{aid,rid});
        ALL_GNDS_EMPTYSPACE{aid, rid} = [x y z xl yl zl];
    end
end


%% representing point
ALL_GNDS_REPPOINT = cell(length(ALL_GNDS), 4);
ALL_GNDS_RULES = cell(length(ALL_GNDS), 4);
for aid = 1:length(ALL_GNDS)
    if isempty(ALL_GNDS{aid})
        continue;
    end
    for rid = 1:1%4
        [ALL_GNDS_REPPOINT{aid, rid} ALL_GNDS_RULES{aid, rid}] = ...
            objectRepPoint( ALL_GNDS_ROT{aid, rid} );
    end
end

end

