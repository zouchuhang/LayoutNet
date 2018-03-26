function compDDSampling( aid )
%COMPDDSAMPLING Summary of this function goes here
%   Detailed explanation goes here
global config;
bufname = config.bufname;

% load([bufname config.objHypoFile3D],'CAN_POOL');
A = load([bufname config.roomModelFile],'orientation','surfacelabel_ori','colorName_rot');
omap = A.orientation;
gc = A.surfacelabel_ori;
cn = A.colorName_rot;
clear A;

% bufName = sprintf([bufname '/' config.objTestImgFeature], aid);
bufName = [bufname '/' config.objTestImgFeature];
% gndBufName = sprintf([bufname '/' config.objGndImgFeature], aid);
gndBufName = [bufname '/' config.objGndImgFeature];
load([bufname config.vpEstimationFile]);
load(config.annotationfile);
H2G_R = ANNO_ALL(aid).ANNO3D.R * ANNO_ALL(aid).ANNO3D.Rc / R;
clear ANNO_ALL;

load(config.groundtruthfile);
gnd = ALL_GNDS_ROT{aid,1};
GND_COMP = getHypInfo({gnd});
clear ALL_GNDS_ROT
load(config.roomtransformfile);
VALID_TRANS_GNDS = TRANS_GNDS(config.anno_valid);
clear TRANS_GNDS

load([bufname config.objHypoFile3D]);
ALL_POOL = CAN_POOL;
clear CAN_POOL

rot_file_name = [bufname config.IMGLIST(aid).name config.rotateImageSuffix];
rotImg = imread(rot_file_name);

%% sample
for rid = 1:length(ALL_POOL)
    fprintf('Room ID: %d\n', rid);
    if config.USEPREVIOUS
        if exist(sprintf([bufname '/' config.sceneFeatFile], rid), 'file') ...
                && exist(sprintf([bufname '/' config.objFeatFile], rid), 'file')
            continue;
        end
    end
    
    CAN_POOL = ALL_POOL(rid);
    if isempty(CAN_POOL{1})
        objfea = [];
        roomfea = [];
        anglesid = [];
        ALL_HYPS = [];
        if config.UPDATE_TO_DISK
            save(sprintf([bufname '/' config.objFeatFile], rid),'objfea','roomfea','anglesid');
            save(sprintf([bufname '/' config.sceneFeatFile], rid), 'ALL_HYPS');
        end
        fprintf('Empty pool, done.\n');
        continue;
    end

    [objfea, roomfea, anglesid] = compObjHypsFeature( rotImg, omap, gc, cn, CAN_POOL{1}, config.typenum, bufName );
    if config.UPDATE_TO_DISK
        save(sprintf([bufname '/' config.objFeatFile], rid),'objfea','roomfea','anglesid');
    end
    
    CAN_POOL{1}.sel_hyps.objfea = objfea;
    CAN_POOL{1}.room.roomfea = roomfea;
    CAN_POOL{1}.sel_hyps.anglesid = anglesid;
    ALL_HYPS = globalSampling( CAN_POOL{1}.sel_hyps, config.samplenum, config.typenum );
    [ sceneImgFea ] = compSceneHypsFeatureA( ALL_HYPS, CAN_POOL{1}.sel_hyps, CAN_POOL{1}.room, omap, gc, cn );
    for i = 1:config.samplenum
        ALL_HYPS(i).sceneImgFea = sceneImgFea(i,:);
    end

    ALL_SCENE = packupScene(ALL_HYPS, CAN_POOL{1});
    [ MINSCORE, MINTRANS] = compRoomMatchScore( ALL_SCENE, VALID_TRANS_GNDS );
    for i = 1:length(ALL_SCENE)
        ALL_HYPS(i).MINSCORE = MINSCORE(i,:);
        ALL_HYPS(i).MINTRANS = MINTRANS(i,:);
    end

    if strcmp(config.lossfunc, 'cost')
        for i = 1:length(ALL_SCENE)
            if isempty(ALL_SCENE{i})
                ALL_HYPS(i).COST = 4.5;
                continue;
            end
            ALL_HYPS(i).COST = roomLossFunction3D( gnd, ALL_SCENE{i}, H2G_R );
        end 
    elseif strcmp(config.lossfunc, 'align')
        for i = 1:length(ALL_SCENE)
            if isempty(ALL_SCENE{i})
                ALL_HYPS(i).COST = 4.5;
                continue;
            end
            HYP_COMP = getHypInfo(ALL_SCENE{i});
            ALL_HYPS(i).COST = roomAlignmentMex(HYP_COMP, GND_COMP, 3, 0);
        end 
    end
    if config.UPDATE_TO_DISK
        save(sprintf([bufname '/' config.sceneFeatFile], rid), 'ALL_HYPS');       
    end
end

% sample near ground truth, only for training
if config.SAMPLE_CLOSE_GND
    [GOOD_HYPS, BESTRHID] = gndtruthSampling( ALL_POOL, gnd, config.samplenum );
    load(sprintf([bufname '/' config.objFeatFile], BESTRHID));
    ALL_POOL{BESTRHID}.sel_hyps.objfea = objfea;
    ALL_POOL{BESTRHID}.room.roomfea = roomfea;
    ALL_POOL{BESTRHID}.sel_hyps.anglesid = anglesid;
    [ sceneImgFea ] = compSceneHypsFeatureA( GOOD_HYPS, ALL_POOL{BESTRHID}.sel_hyps, ALL_POOL{BESTRHID}.room, omap, gc, cn );
    for i = 1:config.samplenum
        GOOD_HYPS(i).sceneImgFea = sceneImgFea(i,:);
    end

    ALL_SCENE = packupScene( GOOD_HYPS, ALL_POOL{BESTRHID});
    [ MINSCORE, MINTRANS] = compRoomMatchScore( ALL_SCENE, VALID_TRANS_GNDS );
    for j = 1:length(ALL_SCENE)
        GOOD_HYPS(j).MINSCORE = MINSCORE(j,:);
        GOOD_HYPS(j).MINTRANS = MINTRANS(j,:);
    end

    if strcmp(config.lossfunc, 'cost')
        for i = 1:length(ALL_SCENE)
            if isempty(ALL_SCENE{i})
                GOOD_HYPS(i).COST = 4.5;
                continue;
            end
            GOOD_HYPS(i).COST = roomLossFunction3D( gnd, ALL_SCENE{i}, H2G_R );
        end 
    elseif strcmp(config.lossfunc, 'align')
        for i = 1:length(ALL_SCENE)
            if isempty(ALL_SCENE{i})
                GOOD_HYPS(i).COST = 4.5;
                continue;
            end
            HYP_COMP = getHypInfo(ALL_SCENE{i});
            GOOD_HYPS(i).COST = roomAlignmentMex(HYP_COMP, GND_COMP, 3, 0);
        end 
    end

    % compute ground truth
    gnd_rotImg = imread([bufname '/' config.IMGLIST(aid).name config.gndImageSuffix]);
    gnd_omap = rotatePanorama( omap, H2G_R);
    gnd_gc = rotatePanorama( gc, H2G_R);
    gnd_cn = rotatePanorama( cn, H2G_R);

    [ GND_HYPS, ~ ] = groundTruthHypothesisNewFunc( gnd, gnd_omap, gnd_gc, gnd_cn, gnd_rotImg, config.typenum, gndBufName);
    [ MINSCORE, MINTRANS] = compRoomMatchScore( {gnd}, VALID_TRANS_GNDS );
    GND_HYPS.MINSCORE = MINSCORE;
    GND_HYPS.MINTRANS = MINTRANS;
    GND_HYPS.COST = 0;
    if config.UPDATE_TO_DISK
        save([bufname config.goodSceneFeatFile], 'GOOD_HYPS', 'GND_HYPS', 'BESTRHID');
    end
end

end

