function [CAN_POOL] = generateObjHyps( aid, numhyp, rect_score )
%GENERATESEEDS Popup object hypothesis from 2D to 3D
%   Detailed explanation goes here
global config


%% load data
rectangle = [];
rectCuboid2 = [];

load([config.bufname config.roomHypoFile]);
load([config.bufname config.objHypoFileUpdate]);
% load([config.folderName config.IMGLIST(aid).name '/' config.roomHypoFile]);
% load([config.folderName config.IMGLIST(aid).name '/' config.objHypoFileUpdate]);
% load([config.folderName config.IMGLIST(aid).name '/' config.rectangleScore]);

load(config.classifierName);
if ~exist('ClassNames','var')
    ClassNames = B.ClassNames;
end
% load('./object_classifier/objectclassifier.mat');
% load('./living_room/objectClassifier_livingroom.mat');

%% update
% ALL_HYPS = cell(numhyp, numsample);

all_opScr = [hypothesis.opScr];
[S,I] = sort(all_opScr, 'descend');
hypothesis = hypothesis(I(1:config.sampleroomnum));
hypothesis_ori = hypothesis_ori(I(1:config.sampleroomnum));

CAN_POOL = cell(min(numhyp, length(hypothesis)),1);

for hid = 1:min(numhyp, length(hypothesis))
% while hid<=numhyp
    room_hyp = hypothesis(hid);
    camera_height = 160;
    [ room_obj ] = room3DFromHyp( room_hyp.hCorner, -camera_height );
    
    room_obj.omapScr = hypothesis(hid).omapScr;
    room_obj.gcScr = hypothesis(hid).gcScr;
    room_obj.mgScr = hypothesis(hid).mgScr;
    room_obj.opScr = hypothesis(hid).opScr;
    
    out_sign = sign(room_obj.out_points_w);
    if ~all(sum(out_sign,1)==0)
        continue;
    end
    
    %% get all hypotheses and project to 3D   
    rectHyps = repmat(struct('objects',[],'scores',[],'valid',[],'bottomscore',[]), 14, 1);

    fprintf('Generaing cuboid: %d/%d, part 1\n', aid, hid);
    parfor rid = 1:5 %parfor
        objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
        bottomscore = zeros(0,1);
        num = size(rectangle(rid).xyzBox,1);
        validIDs = find(rect_score{rid}>0.6);
        for bid = validIDs % 1:num
    %         fprintf('%d/%d\n', bid, num);
            if rid<=4
                newObj = rectangle2hypothesis( rectangle(rid).xyzBox(bid,:), rid, room_obj.out_points_w );
            elseif rid==5
                newObj = floorRectangle2hypothesis( rectangle(rid).xyzBox(bid,:), rid, room_obj.out_points_w );
            end
            objects(end+1:end+length(newObj)) = newObj;
            bottomscore(end+1:end+length(newObj),1) = rectangle(rid).score(bid);
    %         figure(1); clf; viewSingleRectangle(rotImg, rectangle, rid, bid);
    %         figure(2); clf; visualObject3D( [newObj;room_obj], 4, [1 1 1] );
        end
        rectHyps(rid).objects = objects;
        rectHyps(rid).bottomscore = bottomscore;
    end

    fprintf('Generaing cuboid: %d/%d, part 2\n', aid, hid);
    objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
    bottomscore = zeros(0,1);
    for cid = 1:regionCuboid.count
    %     fprintf('%d/%d\n', cid, 100);
        views = regionCuboid.views(cid,:);
        if any(views==6)
            continue;
        end
        rects = regionCuboid.xyzBox(cid,:);
        rects = [rects(1:12); rects(13:24); rects(25:36)];
    %     newObj = obj3DFromHyp( rects, views, room_obj.out_points_w ); 
        newObj = rectangleCombo3Hypothesis( rects, views, room_obj.out_points_w );
        objects(end+1:end+length(newObj)) = newObj;
        bottomscore(end+1:end+length(newObj),1) = regionCuboid.score(cid);
    %     figure(1); clf; viewRegionHypothesisSubset( rotImg, regionCuboid, cid );
    %     figure(2); clf; visualObject3D( [newObj;room_obj], 4 );   
    end
    rectHyps(6).objects = objects;
    rectHyps(6).bottomscore = bottomscore;

    fprintf('Generaing cuboid: %d/%d, part 3\n', aid, hid);
    parfor rid = 1:8 %parfor
    %     fprintf('%d/%d\n', rid, length(rectCuboid2));
        views = rectCuboid2(rid).views;
        views(views==0) = [];
        score1 = rectangle(views(1)).score;
        score2 = rectangle(views(2)).score;
        if any(views==5)
            num = length(rectCuboid2(rid).score);
            objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
            bottomscore = zeros(0,1);
            for cid = 1:num
                combo = rectCuboid2(rid).combo(cid,:);
                rects = [rectangle(views(1)).xyzBox(combo(1),1:12); rectangle(views(2)).xyzBox(combo(2),1:12)];
                newObj = floorRectangleCombo2Hypothesis( rects, views, room_obj.out_points_w );

    %             figure(1); clf; viewRectHypothesisSubset(rotImg, rectCuboid2, rectangle, rid, cid);
    %             figure(2); clf; visualObject3D( [newObj;room_obj], 4 );
                objects(end+1:end+length(newObj)) = newObj;
                bottomscore(end+1:end+length(newObj),1) = max(score1(cid), score2(cid)); %rectCuboid2(rid).score(cid);
            end
        else
            num = length(rectCuboid2(rid).score);
            objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
            bottomscore = zeros(0,1);
            for cid = 1:num
                combo = rectCuboid2(rid).combo(cid,:);
                rects = [rectangle(views(1)).xyzBox(combo(1),1:12); rectangle(views(2)).xyzBox(combo(2),1:12)];
                newObj = rectangleCombo2Hypothesis( rects, views, room_obj.out_points_w );

    %             figure(1); clf; viewRectHypothesisSubset(rotImg, rectCuboid2, rectangle, rid, cid);
    %             figure(2); clf; visualObject3D( [newObj;room_obj], 4 );
                objects(end+1:end+length(newObj)) = newObj;
                bottomscore(end+1:end+length(newObj),1) = max(score1(cid), score2(cid)); %rectCuboid2(rid).score(cid);
            end
        end

        rectHyps(rid+6).objects = objects;
        rectHyps(rid+6).bottomscore = bottomscore;
    end
    
    %% classify: semantic classifier, should move to DDsampling   
    for rid = 1:14
        [ scores, valid ] = objectClassifier(room_obj, rectHyps(rid).objects, B, ClassNames);
        rectHyps(rid).scores = scores;
        rectHyps(rid).valid = valid;
%         num = length(rectHyps(rid).objects);
%         rectHyps(rid).scores = zeros(num,28);
%         rectHyps(rid).valid = false(num,1);
    end
    all_objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
    scores = zeros(0,28);
    valid = false(0,1);
    obj_xyz = zeros(0,6);
    bottomscore = zeros(0,1);
    for rid = 1:14
        num = length(rectHyps(rid).objects);
        all_objects(end+1:end+num,:) = rectHyps(rid).objects;
        scores(end+1:end+num,:) = rectHyps(rid).scores ;
        valid(end+1:end+num,1) = rectHyps(rid).valid;
        bottomscore(end+1:end+num,1) = rectHyps(rid).bottomscore;
        xyz = zeros(num,6);
        for i = 1:num
            xyz(i,:) = [min(rectHyps(rid).objects(i).out_points_w, [], 1) ...
                max(rectHyps(rid).objects(i).out_points_w, [], 1)];
        end
        obj_xyz(end+1:end+num,:) = xyz;
    end
    
    scores(~valid,:) = 0;
    all_hyps.objects = addBoundingBox(all_objects(valid),'both'); %all_objects(valid);
    all_hyps.scores = scores(valid,:);
    all_hyps.valid = valid(valid);
    all_hyps.obj_xyz = obj_xyz(valid,:);
    all_hyps.bottomscore = bottomscore(valid);
    fprintf('%d:%d, %d-->%d\n', aid, hid, length(valid), sum(valid));
 
    
    
    %% prune weak cuboid hypothesis
%     t = 0;
%     sel_hyps = [];
%     for oid = 1:28
%         valid = bottomscore>0 & scores(:,oid)>t;
%         sel_hyps(oid).h_objects = all_hyps.objects(valid);
%         sel_hyps(oid).h_scores = all_hyps.scores(valid,:);
%         sel_hyps(oid).h_obj_xyz = all_hyps.obj_xyz(valid,:);   
%         sel_hyps(oid).h_bottomscore = all_hyps.bottomscore(valid,:);
%         sel_hyps(oid).h_objID = find(valid);
%     end
% %     t = 0.60;
% %     sel_hyps = [];
% %     for oid = 1:28
% %         [S,I] = sort(all_hyps.scores(:,oid),'descend');
% %         obj_xyz = all_hyps.obj_xyz(I(S>t),:);
% %         [ outxyz, valid ] = nms3D( obj_xyz, 0.5 );
% %         sel_hyps(oid).objects = all_hyps.objects(I(valid));
% %         sel_hyps(oid).scores = all_hyps.scores(I(valid),:);
% %         sel_hyps(oid).obj_xyz = all_hyps.obj_xyz(I(valid),:);   
% %         sel_hyps(oid).bottomscore = all_hyps.bottomscore(I(valid),:);
% %         sel_hyps(oid).objID = I(valid);
% % 
% %         v = sel_hyps(oid).bottomscore>0;
% % 
% %         sel_hyps(oid).h_objects = sel_hyps(oid).objects(v);
% %         sel_hyps(oid).h_objID = sel_hyps(oid).objID(v);
% %         sel_hyps(oid).h_scores = sel_hyps(oid).scores(v,:);
% %         sel_hyps(oid).h_obj_xyz = sel_hyps(oid).obj_xyz(v,:);
% %         sel_hyps(oid).h_bottomscore = sel_hyps(oid).bottomscore(v,:);
% %     end
% 
%     % pack gnd 
%     [ALL_HYPS(hid,:), CAN_POOL{hid}.sel_hyps, CAN_POOL{hid}.room] = getInitialWholeLayoutB( sel_hyps, room_obj, numsample );
    CAN_POOL{hid}.sel_hyps = all_hyps;
    CAN_POOL{hid}.room = room_obj;
end

% ALL_HYPS = ALL_HYPS(:);

% %% output candidate pool
% for tid = 1:12
%     
% end

end

