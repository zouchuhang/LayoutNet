function [ output_args ] = panoContext( aid )
%PANOCONTEXT Summary of this function goes here
%   Detailed explanation goes here
global config;
UPDATE_TO_DISK = false;
bufname = [config.folderName config.IMGLIST(aid).name '/'];

fprintf('Annotation ID: %d\n', aid);

%% vanishing point
COMPUTE_VP = false;

if COMPUTE_VP
    Img = imread(config.IMGLIST(aid).imgname);
    Img = im2double(Img);
    imgSize = 320;
    qError = 0.7;
    [ olines, vp, views, edges, panoEdge, score, angle] ...
        = panoEdgeDetection( Img, imgSize, qError); 

    Img_small = imresize(Img, [1024 2048]);
    [rotImg, R] = rotatePanorama(Img_small, vp(3:-1:1,:));

    if UPDATE_TO_DISK
        sml_file_name = [bufname config.IMGLIST(aid).name config.smallImageSuffix];
        imwrite(Img_small, sml_file_name);

        rot_file_name = [bufname config.IMGLIST(aid).name config.rotateImageSuffix];
        imwrite(rotImg, rot_file_name);

        parsave([bufname config.vpEstimationFile], 'vp', vp, 'R', R);
    end
else
    sml_file_name = [bufname config.IMGLIST(aid).name config.smallImageSuffix];
    smlImg = im2double(imread(sml_file_name));
    rot_file_name = [bufname config.IMGLIST(aid).name config.rotateImageSuffix];
    rotImg = im2double(imread(rot_file_name));
    load([bufname config.vpEstimationFile]);
end
%% sample room hypothesis
COMPUTE_ROOM = false;

if COMPUTE_ROOM
    [ ~, panoOmap ] = computePanoOmap( views, edges, vp );
    panoOmap_rot = rotatePanorama(panoOmap, [], R);   
    [ ~, wallPanoNormal] = compSurfaceLabel( rotImg );
    wallPanoNormal_rot = rotatePanorama(wallPanoNormal, [], R);

    orientation_ori = panoOmap; %C.orientation_ori;
    surfacelabel = rotatePanorama(wallPanoNormal, inv(R)); %rotatePanorama(C.surfacelabel_ori, inv(R));
    hyps = generateHypsB(olines, vp, 3, orientation_ori, surfacelabel);
    hyps_rot = rotateHyps( hyps, R);
    
    colorName_rot = im2cM(double(rotImg));
    colorName = rotatePanorama(colorName_rot, inv(R));

    if UPDATE_TO_DISK
        parsave([bufname config.roomModelFile], ...
            'orientation_ori', panoOmap, 'surfacelabel_ori', wallPanoNormal, ...
            'orientation', panoOmap_rot, 'surfacelabel', wallPanoNormal_rot, ...
            'colorName', colorName, 'colorName_rot', colorName_rot);
        parsave([bufname config.roomHypoFile], ...
            'hypothesis_ori', hyps, 'hypothesis', hyps_rot);
    end
else
    load([bufname config.roomHypoFile]);
    load([bufname config.roomModelFile]);
end
    
%% generate object hypothesis
COMPUTE_OBJECT = false;

if COMPUTE_OBJECT
    % rectangle hypothesis
    [ rectangle, rectCuboid2, rectCuboid3 ] = rectangleBasedHypothesis( rotImg );   
    % segmentation hypothesis
    [ regionCuboid ] = regionBasedHypothesis( rotImg );
    % selective search
    SS = getSelectiveSearch(rotImg);
    AllCuboid = regionBasedHypothesisFromSS(SS, 1024, 2048);
    rectangle = rectDetectionFlexFromSS(SS, 1024, 2048);
    AllCuboid_NEW = pruneCuboid(rotImg, AllCuboid);
    rectangle_NEW = pruneRectangle(rotImg, rectangle);
    % merge up everything
    regionCuboid.views = [regionCuboid.views; AllCuboid_NEW.views];
    regionCuboid.xyzBox = [regionCuboid.xyzBox; AllCuboid_NEW.xyzBox];
    regionCuboid.score = [regionCuboid.score; AllCuboid_NEW.score];
    regionCuboid.count = regionCuboid.count + AllCuboid_NEW.count;

    for i = 1:6
        rectangle(i).highPixelBox = [rectangle(i).highPixelBox; rectangle_NEW(i).highPixelBox];
        rectangle(i).xyzBox = [rectangle(i).xyzBox; rectangle_NEW(i).xyzBox];
        rectangle(i).score = [rectangle(i).score; rectangle_NEW(i).score];
    end

    if UPDATE_TO_DISK
        parsave([bufname config.objHypoFileUpdate], ...
            'rectangle', rectangle, 'rectCuboid2', rectCuboid2, 'rectCuboid3', rectCuboid3, ...
            'regionCuboid', regionCuboid);
    end

    [ rect_score ] = getSegScrOfRect( rotImg, rectangle );
    save([bufname config.rectangleScore], 'rect_score');
    CAN_POOL = generateObjHyps( aid, 20000);
    save([bufname config.objHypoFile3D],'CAN_POOL', '-v7.3');
    fprintf('finish.\n');
end

%% data driven sampling
COMPUTE_SAMPLING = false;
ALL_POOL = CAN_POOL;
if COMPUTE_SAMPLING
    omap = panoOmap_rot;
    gc = wallPanoNormal;
    cn = colorName_rot;
    bufName = sprintf([bufname '/' config.objTestImgFeature], aid);
    gndBufName = sprintf([bufname '/' config.objGndImgFeature], aid);
    H2G_R = ANNO_ALL(aid).ANNO3D.R * ANNO_ALL(aid).ANNO3D.Rc / R;
    % sample on image
    for rid = 1:length(ALL_POOL)
        CAN_POOL = ALL_POOL(rid);
        if isempty(POOL{1})
            objfea = [];
            roomfea = [];
            anglesid = [];
            ALL_HYPS = [];
            if UPDATE_TO_DISK
                save(sprintf([bufname '/' config.objFeatFile], rid),'objfea','roomfea','anglesid');
                save(sprintf([bufname '/' config.sceneFeatFile], rid), 'ALL_HYPS');
            end
            fprintf('Empty pool, done.\n');
            continue;
        end
       
        [objfea, roomfea, anglesid] = compObjHypsFeature( rotImg, omap, gc, cn, CAN_POOL{1}, config.typenum, bufName );
        save(sprintf([bufname '/' config.objFeatFile], rid),'objfea','roomfea','anglesid');
        
        CAN_POOL{1}.sel_hyps.objfea = objfea;
        CAN_POOL{1}.room.roomfea = roomfea;
        CAN_POOL{1}.sel_hyps.anglesid = anglesid;
        ALL_HYPS = globalSampling( CAN_POOL{1}.sel_hyps, CAN_POOL{1}.room, config.samplenum, config.typenum );
        [ sceneImgFea ] = compSceneHypsFeature( ALL_HYPS, CAN_POOL{1}.sel_hyps, CAN_POOL{1}.room, omap, gc, cn );
        for i = 1:num_sample
            ALL_HYPS(i).sceneImgFea = sceneImgFea(i,:);
        end
        
        ALL_SCENE = packupScene(ALL_HYPS, CAN_POOL{1});
        [ MINSCORE, MINTRANS] = compRoomMatchScore( ALL_SCENE, config.VALID_TRANS_GNDS );
        for i = 1:length(ALL_SCENE)
            ALL_HYPS(i,hid).MINSCORE = MINSCORE(i,:);
            ALL_HYPS(i,hid).MINTRANS = MINTRANS(i,:);
        end
       
        for i = 1:length(ALL_SCENE)
            if isempty(ALL_SCENE{i})
                ALL_HYPS(i,hid).COST = 4.5;
                continue;
            end
            ALL_HYPS(i,hid).COST = roomLossFunction3D( config.ALL_GNDS_ROT{aid,1}, ALL_SCENE{i}, H2G_R );
        end 
        if UPDATE_TO_DISK
            save(sprintf([bufname '/' config.sceneFeatFile], rid), 'ALL_HYPS');       
        end
    end
    % sample near ground truth, only for training
    [GOOD_HYPS, BESTRHID] = gndtruthSampling( ALL_POOL, config.ALL_GNDS_ROT{aid,1}, config.samplenum );
    load(sprintf([bufname '/' config.objFeatFile], BESTRHID));
    ALL_POOL{BESTRHID}.sel_hyps.objfea = objfea;
    ALL_POOL{BESTRHID}.room.roomfea = roomfea;
    ALL_POOL{BESTRHID}.sel_hyps.anglesid = anglesid;
    [ sceneImgFea ] = compSceneHypsFeature( GOOD_HYPS, ALL_POOL{BESTRHID}.sel_hyps, ALL_POOL{BESTRHID}.room, omap, gc, cn );
    for i = 1:num_sample
        GOOD_HYPS(i).sceneImgFea = sceneImgFea(i,:);
    end

    ALL_SCENE = packupScene( GOOD_HYPS, ALL_POOL{BESTRHID});
    [ MINSCORE, MINTRANS] = compRoomMatchScore( ALL_SCENE, config.VALID_TRANS_GNDS );
    for j = 1:length(ALL_SCENE)
        GOOD_HYPS(j).MINSCORE = MINSCORE(j,:);
        GOOD_HYPS(j).MINTRANS = MINTRANS(j,:);
    end

    for j = 1:length(ALL_SCENE)
        if isempty(ALL_SCENE{j})
            GOOD_HYPS(j).COST = 4.5;
            continue;
        end
        GOOD_HYPS(j).COST = roomLossFunction3D( config.ALL_GNDS_ROT{aid,1}, ALL_SCENE{j}, H2G_R );
    end 
    % compute ground truth
    gnd_rotImg = imread([bufname '/' config.IMGLIST(aid).name config.gndImageSuffix]);
    gnd_omap = rotatePanorama( omap, H2G_R);
    gnd_gc = rotatePanorama( gc, H2G_R);
    gnd_cn = rotatePanorama( cn, H2G_R);
    
    [ GND_HYPS, ~ ] = groundTruthHypothesisNewFunc( config.ALL_GNDS_ROT{aid,1}, gnd_omap, gnd_gc, gnd_cn, gnd_rotImg, config.typenum, gndBufName);
    [ MINSCORE, MINTRANS] = compRoomMatchScore( config.ALL_GNDS_ROT(aid,1), config.VALID_TRANS_GNDS );
    GND_HYPS.MINSCORE = MINSCORE;
    GND_HYPS.MINTRANS = MINTRANS;
    GND_HYPS.COST = 0;
    if UPDATE_TO_DISK
        save([folderName ANNO_ALL(aid).name '/' config.goodSceneFeatFile], 'GOOD_HYPS', 'GND_HYPS', 'BESTRHID');
    end
end

%% holistic ranking



end

