function [ ALL_HYPS, BESTRHID ] = gndtruthSampling( CAN_POOL, GT, number )
%GNDTRUTHSAMPLING Summary of this function goes here
%   Detailed explanation goes here

% v = ~cellfun(@isempty, CAN_POOL);
% CAN_POOL = CAN_POOL(v);
for i = 1:length(CAN_POOL)
    if isempty(CAN_POOL{i})
        continue;
    end
    CAN_POOL{i}.room_E = addBoundingBox(CAN_POOL{i}.room, 'both');  
end

%% find closest room hypothesis
roomlo_score = zeros(length(CAN_POOL),1);
for i = 1:length(CAN_POOL)
    if isempty(CAN_POOL{i})
        roomlo_score(i) = 10000000;
        continue;
    end
    pt_gt = GT(1).align;
    pt_hp = CAN_POOL{i}.room_E.align;
    roomlo_score(i) = sum((pt_gt(:)-pt_hp(:)).^2);
end
[~,BESTRHID] = min(roomlo_score);

%% sample closest room
room_obj = CAN_POOL{BESTRHID}.room;
type_hyps = CAN_POOL{BESTRHID}.sel_hyps;
typenum = 12;
for tid = 1:typenum
    sel_hyps(tid).fixed = false;
    sel_hyps(tid).selID = [];
end
goodIDs = cell(length(GT),1);

for oid = 2:length(GT)
    objtype = GT(oid).objtype;
    if objtype>12
        continue;
    end
    hypnum = length(type_hyps.objects);
    if hypnum==0
        continue;
    end
    objscore = zeros(hypnum, 1);
    for hid = 1:hypnum
        objscore(hid) = sum((type_hyps.objects(hid).align(:)-GT(oid).align(:)).^2);
    end
    [s,bestobjid] = sort(objscore);
    bestobjid = bestobjid(s<10^5);
    goodIDs{oid} = bestobjid;
    
    count = 1;
    while count<=length(bestobjid) && any(sel_hyps(objtype).selID==bestobjid(count))
        count = count + 1;
    end
    
    if count>length(bestobjid)
        continue;
    end
    
    if ~isempty(bestobjid)
        sel_hyps(objtype).fixed = true;
        sel_hyps(objtype).selID = [sel_hyps(objtype).selID bestobjid(count)];
    end
end

ALL_HYPS(1).sel_hyps = sel_hyps;
% [best_objects] = visualWholeRoomObject3D( type_hyps, sel_hyps, 1, room_obj);

%% sample closer room

for sid = 2:number
    for tid = 1:typenum
        sel_hyps(tid).fixed = false;
        sel_hyps(tid).selID = [];
    end
    for oid = 2:length(GT)
        objtype = GT(oid).objtype;  
        if isempty(goodIDs{oid})
            continue;
        end
        bestobjid = randi(min(10, length(goodIDs{oid})),1);
        randid = goodIDs{oid}(bestobjid);
        sel_hyps(objtype).fixed = true;
        count=0;
        while any(sel_hyps(objtype).selID==randid) && count<10
            bestobjid = randi(min(10, length(goodIDs{oid})),1);
            randid = goodIDs{oid}(bestobjid);
            count = count + 1;
        end
        if count==10
            continue;
        end
        sel_hyps(objtype).selID = [sel_hyps(objtype).selID randid];
    end
%     [objects] = visualWholeRoomObject3D( type_hyps, sel_hyps, 1, room_obj);
    ALL_HYPS(sid).sel_hyps = sel_hyps;
end
% cost = zeros(number, 1);
% for i = 1:number
%     cost(i) = roomLossFunction3D( GT, ALL_HYPS{i}, diag([1 1 1]) );
% end

end

