function compLocalSampling( aid )
%COMPHOLISTICRANKING Summary of this function goes here
%   Detailed explanation goes here
global config
bufname = config.bufname;

did = config.valid_data_id(config.valid_anno_id==aid);

load(config.globalSceneSVMFile);
scenesvm = MODEL{config.globalSceneSVMiterID}.initModel;
clear MODEL

load(sprintf([bufname config.sceneFeatFile], 1));
SEED_HYPS = ALL_HYPS(1);
names = fieldnames(SEED_HYPS);
for i = 1:length(names)
    SEED_HYPS.(names{i}) = [];
end

HYPS = repmat( SEED_HYPS, 5000, 1);
COST = 4.5*ones(5000,1);
ROOMID = zeros(5000,1);
SCORE = -10*ones(5000,1);
for rid = 1:config.sampleroomnum
    fprintf('>>%d:%d ', aid, rid);
    if exist(sprintf([bufname config.sceneFeatFile], rid), 'file')
        load(sprintf([bufname config.sceneFeatFile], rid));
    else
        continue;
    end
    if isempty(ALL_HYPS)
        continue;
    end
    
    HYPS((rid-1)*100+1:rid*100) = ALL_HYPS;
    if strcmp(config.lossfunc,'align') && isfield(ALL_HYPS.COSTALIGN) % cope with current version, later there would be no COSTALIGN
        COST((rid-1)*100+1:rid*100) = vertcat(ALL_HYPS.COSTALIGN);
    else
        COST((rid-1)*100+1:rid*100) = vertcat(ALL_HYPS.COST);
    end
    ROOMID((rid-1)*100+1:rid*100) = rid;
    
    [ ~, hypScore] = sceneHypsEvaluation( vertcat(ALL_HYPS.sceneImgFea), ...
        vertcat(ALL_HYPS.MINSCORE), scenesvm, did );
    SCORE((rid-1)*100+1:rid*100) = hypScore;  
end

[~,I] = sort(SCORE,'descend');

res.BESTHYPS = HYPS(I(1:50));
res.BESTROOM = ROOMID(I(1:50));
res.BESTSCRS = SCORE(I(1:50));
res.BESTCOST = COST(I(1:50));
save([bufname config.globalResultFile], 'res');

%% local sampling
A = load([bufname config.roomModelFile],'orientation','surfacelabel_ori','colorName_rot');
omap = A.orientation;
gc = A.surfacelabel_ori;
cn = A.colorName_rot;
clear A;
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


ROOM_IDS = res.BESTROOM(1:3);
ROOM_IDS = unique(ROOM_IDS);
ITERHYPS = cell(length(ROOM_IDS),10);

load([bufname config.objHypoFile3D]);
load(config.numberName);
for uid = 1:length(ROOM_IDS)
    roomid = ROOM_IDS(uid);
    fprintf('ROOM ID: %d\n', roomid);
    load(sprintf([bufname config.objFeatFile], roomid)); 
    
    OBJ_POOL = CAN_POOL(roomid);
    
    try
        unaryscore = vertcat(objfea.unaryscore);
%         unaryscore = additionalGeoRule(OBJ_POOL{1}.sel_hyps.obj_xyz, unaryscore, 1000000);
        unaryscore = feval(config.georulefunc, OBJ_POOL{1}.sel_hyps.obj_xyz, unaryscore, 1000000);
        for i = 1:length(objfea)
            objfea(i).unaryscore = unaryscore(i,:);
        end
    catch
        fprintf('%d, %d\n', aid, roomid);
        fprintf('%d, %d', size(unaryscore,1), length(objfea));
        fprintf('%d', anything);
    end
      
    OBJ_POOL{1}.sel_hyps.objfea = objfea;
    OBJ_POOL{1}.sel_hyps.anglesid = anglesid;
    OBJ_POOL{1}.room.roomfea = roomfea;
    
    load(sprintf([bufname config.sceneFeatFile], roomid));
       
    obj_xyz = OBJ_POOL{1}.sel_hyps.obj_xyz;
    % evaluate all objects
    [ objScore ] = quickObjectEvalA( objfea, scenesvm, config.typenum );
    objScore = 1./(1+exp(-0.5*(objScore+6)));
    numobj = length(objScore);
    NEW_HYPS = ALL_HYPS;
    WHOLE_HYPS = repmat(NEW_HYPS(1), 0, 1);
    WHOLE_SCRS = zeros(0,1); 
    
    for iter = 1:100
        fprintf('DATA: %d, UID: %d, ITERATION: %d\n', aid, uid, iter);        
        fprintf('Evaluate current hypotheses set\n');
        % evaluate all scenes
        [ ~, hypScore] = sceneHypsEvaluation( vertcat(NEW_HYPS.sceneImgFea), ...
                vertcat(NEW_HYPS.MINSCORE), scenesvm, did );
        WHOLE_HYPS = [WHOLE_HYPS; NEW_HYPS];
        WHOLE_SCRS = [WHOLE_SCRS; hypScore];
            
        [~, I] = sort(WHOLE_SCRS, 'ascend');
        J = 1:length(I); J(I) = J; J = J' - (length(I)-25);
        selvalid = softSelect( J, 0, 100 );        
        I = find(selvalid); 
        [~,J] = sort(WHOLE_SCRS(I),'descend');
        hypIDs = I(J);
        SEL_HYPS = WHOLE_HYPS(hypIDs);
        
        ITERHYPS{uid,iter}.Hyps = SEL_HYPS;
        ITERHYPS{uid,iter}.Scrs = WHOLE_SCRS(hypIDs);
        ITERHYPS{uid,iter}.roomid = roomid;
        if  iter>=4
            if iter == config.MAX_ITER
                break;
            elseif ITERHYPS{uid,iter}.Scrs(iter)-ITERHYPS{uid,iter-3}.Scrs(1)<=0
                break;
            end
        end
        
        SEL_HYPS = rmfield(SEL_HYPS, {'sceneImgFea','MINSCORE','MINTRANS','COST'});
        
        objvalid = false(length(objScore),config.typenum);  
        for i = 1:length(hypIDs)
            sel_hyps = SEL_HYPS(i).sel_hyps;
            for k = 1:config.typenum
                objvalid(sel_hyps(k).selID,k) = true;
            end
        end  
              
       %% start sampling, top 5, thoroughly sampling
        NEW_HYPS = repmat(SEL_HYPS(1),0,1);
        fprintf('Dense sampling.\n');
        for i = 1:5
            fprintf('*');
            % delete
            sel_hyps = SEL_HYPS(i).sel_hyps;
            num_obj = length([sel_hyps.selID]);
            [ DEL_HYPS, DEL_TYPE, DEL_OBID ] = deleteObjectPerm( 1:num_obj, repmat(SEL_HYPS(i), num_obj, 1) );
            
            % add
            typeprod = zeros(1,config.typenum);
            for j = 1:config.typenum
                typeprod(j) = prob(j,length(sel_hyps(j).selID)+2);%/(prob(j,length(sel_hyps(j).selID)+1)+0.0001);
            end
            TPS = randp(typeprod, 20, 1);
            [ ADD_HYPS_CUR ] = addObject( OBJ_POOL{1}.sel_hyps, repmat(SEL_HYPS(i),20, 1), TPS, ...
                objvalid(:,TPS), objScore ); 
            [ ADD_HYPS_ALL ] = addObject( OBJ_POOL{1}.sel_hyps, repmat(SEL_HYPS(i),20, 1), TPS, ...
                true(numobj,20), objScore ); 
            v = deleteReplicateHyps(ADD_HYPS_CUR);
            ADD_HYPS_CUR = ADD_HYPS_CUR(v);
            v = deleteReplicateHyps(ADD_HYPS_ALL);
            ADD_HYPS_ALL = ADD_HYPS_ALL(v);
                      
            % replace
%             repvalid = false(numobj, 100);
            REP_HYPS_CUR = repmat(SEL_HYPS(1),0,1);
            REP_HYPS_ALL = repmat(SEL_HYPS(1),0,1);
            for j = 1:length(DEL_OBID)
                [ score ] = findNearbyObject( obj_xyz, obj_xyz(DEL_OBID(j),:) );
                repvalid = score>0.2 & score<0.8;
                
                [ rep1 ] = addObject( OBJ_POOL{1}.sel_hyps, repmat(DEL_HYPS(j), 5, 1), ...
                    DEL_TYPE(j)*ones(5,1), ...
                    repmat(repvalid.*objvalid(:,DEL_TYPE(j)),1,5), objScore );
                
                [ rep2 ] = addObject( OBJ_POOL{1}.sel_hyps, repmat(DEL_HYPS(j), 5, 1), ...
                    DEL_TYPE(j)*ones(5,1), ...
                    repmat(repvalid,1,5), objScore );
                
                REP_HYPS_CUR = [REP_HYPS_CUR;rep1];
                REP_HYPS_ALL = [REP_HYPS_ALL;rep2];
            end
            v = deleteReplicateHyps(REP_HYPS_CUR);
            REP_HYPS_CUR = REP_HYPS_CUR(v);
            v = deleteReplicateHyps(REP_HYPS_ALL);
            REP_HYPS_ALL = REP_HYPS_ALL(v);
            
            % delete and add
            DEL_ADD_CUR = repmat(SEL_HYPS(1),0,1);
            DEL_ADD_ALL = repmat(SEL_HYPS(1),0,1);
%             for j = 1:length(DEL_OBID)
%                 hyps = DEL_HYPS(j).sel_hyps;
%                 typeprod = zeros(1,12);
%                 for k = 1:12
%                     typeprod(k) = prob(k,length(hyps(k).selID)+2);%/(prob(k,length(hyps(k).selID)+1)+0.0001);
%                 end
%                 TPS = randp(typeprod, 20, 1);
%                 [ add1 ] = addObject( OBJ_POOL{1}.sel_hyps, repmat(DEL_HYPS(j),20, 1), TPS, ...
%                     objvalid(:,TPS), objScore ); 
%                 [ add2 ] = addObject( OBJ_POOL{1}.sel_hyps, repmat(DEL_HYPS(j),20, 1), TPS, ...
%                     true(numobj,20), objScore ); 
%                 DEL_ADD_CUR = [DEL_ADD_CUR; add1];
%                 DEL_ADD_ALL = [DEL_ADD_ALL; add2];
%             end
            
            NEW_HYPS = [NEW_HYPS; DEL_HYPS; ADD_HYPS_CUR; ADD_HYPS_ALL; ...
                REP_HYPS_CUR; REP_HYPS_ALL; DEL_ADD_CUR; DEL_ADD_ALL];
        end
        v = deleteReplicateHyps([SEL_HYPS(1:5);NEW_HYPS]);
        NEW_HYPS = NEW_HYPS(v(6:end));
        fprintf('\n');
        
       %% for the next, random 1
        fprintf('Sparse sampling\n');
        ALL_HYPS = SEL_HYPS(6:min(length(SEL_HYPS), 50));
        v = true(length(ALL_HYPS),1);
        for i = 1:length(ALL_HYPS)
            o = [ALL_HYPS(i).sel_hyps.selID];
            if isempty(o)
                v(i) = false;
            end
        end
        ALL_HYPS = ALL_HYPS(v);
        
        % step 1: randomly delete an object 
        fprintf('Random move by deleting object.\n');
        [ DEL_HYPS, ~, ~ ] = deleteObject( objScore, ALL_HYPS );
        
        % step 2: randomly add an object with current pool
        fprintf('Random move by adding object from current pool.\n');
        TPS = zeros(length(ALL_HYPS),1);
        for i = 1:length(ALL_HYPS)
            sel_hyps = ALL_HYPS(i).sel_hyps;
            typeprob = zeros(1,config.typenum);
            for j = 1:config.typenum
                typeprob(j) = prob(j,length(sel_hyps(j).selID)+2)/(prob(j,length(sel_hyps(j).selID)+1)+0.0001);
            end
            TPS(i) = randp(typeprob, 1);
        end
        
%         randnum = randp(objScore(:), length(ALL_HYPS), 1);
%         [~, TPS] = ind2sub([numobj typenum], randnum);
        [ ADD_HYPS_CUR ] = addObject( OBJ_POOL{1}.sel_hyps, ALL_HYPS, TPS, ...
            objvalid(:,TPS), objScore ); 
        
        % step 2.5: randomly add an object with all pool
        fprintf('Random move by adding object from all pool.\n');
%         randnum = randp(objScore(:), length(ALL_HYPS),1);
%         [~, TPS] = ind2sub([numobj typenum], randnum);
        [ ADD_HYPS_ALL ] = addObject( OBJ_POOL{1}.sel_hyps, ALL_HYPS, TPS, ...
            true(numobj,length(ALL_HYPS)), objScore ); 
        
        % step 3: randomly replace an object with current pool
        fprintf('Random move by replacing object from current pool.\n');
        [ DEL_HYPS_REP, DEL_TYPE, DEL_OBID ] = deleteObject( objScore, ALL_HYPS );
        repvalid = false(numobj, length(ALL_HYPS));
        for i = 1:length(ALL_HYPS)
            [ score ] = findNearbyObject( obj_xyz, obj_xyz(DEL_OBID(i),:) );
            repvalid(:,i) = score>0.2 & score<0.8;
        end
        [ REP_HYPS_CUR ] = addObject( OBJ_POOL{1}.sel_hyps, DEL_HYPS_REP, DEL_TYPE, ...
            repvalid.*objvalid(:,DEL_TYPE), objScore ); 
        
        % step 3.5: randomly replace an object with all pool
        fprintf('Random move by replacing object from all pool.\n');
        [ DEL_HYPS_REP, DEL_TYPE, DEL_OBID ] = deleteObject( objScore, ALL_HYPS );
        repvalid = false(numobj, length(ALL_HYPS));
        for i = 1:length(ALL_HYPS)
            [ score ] = findNearbyObject( obj_xyz, obj_xyz(DEL_OBID(i),:) );
            repvalid(:,i) = score>0.2 & score<0.8;
        end
        [ REP_HYPS_ALL ] = addObject( OBJ_POOL{1}.sel_hyps, DEL_HYPS_REP, DEL_TYPE, ...
            repvalid, objScore ); 
        
        % step 4: randomly sample with current pool
        fprintf('Random sample from current pool.\n');
        [ RND_HYPS_CUR ] = globalSamplingWhitelist( OBJ_POOL{1}.sel_hyps, ...
            OBJ_POOL{1}.room, 20, config.typenum, objvalid, objScore );  
        
        % step 4.5: randomly sample with all
        fprintf('Random sample from all pool.\n');
        [ RND_HYPS_ALL ] = globalSamplingWhitelist( OBJ_POOL{1}.sel_hyps, ...
            OBJ_POOL{1}.room, 20, config.typenum, true(numobj,config.typenum), objScore ); 
        
        NEW_HYPS = [NEW_HYPS; DEL_HYPS; ADD_HYPS_CUR; ADD_HYPS_ALL; ...
            REP_HYPS_CUR; REP_HYPS_ALL; RND_HYPS_CUR; RND_HYPS_ALL];
        
        NUM = min(length(SEL_HYPS), 50);
        v = deleteReplicateHyps([SEL_HYPS(1:NUM);NEW_HYPS]);
        NEW_HYPS = NEW_HYPS(v(NUM+1:end));
        
        fprintf('Computing scene feature for %d newly sampled hypotheses.\n', length(NEW_HYPS));
        [ sceneImgFea ] = compSceneHypsFeatureA( NEW_HYPS, ...
            OBJ_POOL{1}.sel_hyps, OBJ_POOL{1}.room, omap, gc, cn );
        for i = 1:length(NEW_HYPS)
            NEW_HYPS(i).sceneImgFea = sceneImgFea(i,:);
        end
        fprintf('Computing room matching cost for %d newly sampled hypotheses.\n', length(NEW_HYPS));
        ALL_SCENE = packupScene(NEW_HYPS, OBJ_POOL{1});
        [ MINSCORE, MINTRANS] = compRoomMatchScore( ALL_SCENE, VALID_TRANS_GNDS );
        for i = 1:length(ALL_SCENE)
            NEW_HYPS(i).MINSCORE = MINSCORE(i,:);
            NEW_HYPS(i).MINTRANS = MINTRANS(i,:);           
        end
        fprintf('Computing room gnd cost for %d newly sampled hypotheses.\n', length(NEW_HYPS));
        if strcmp(config.lossfunc, 'cost')
            for i = 1:length(ALL_SCENE)
                if isempty(ALL_SCENE{i})
                    NEW_HYPS(i).COST = 4.5;
                    continue;
                end
                NEW_HYPS(i).COST = roomLossFunction3D( gnd, ALL_SCENE{i}, H2G_R );
            end 
        elseif strcmp(config.lossfunc, 'align')
            for i = 1:length(ALL_SCENE)
                if isempty(ALL_SCENE{i})
                    NEW_HYPS(i).COST = 4.5;
                    continue;
                end
                HYP_COMP = getHypInfo(ALL_SCENE{i});
                NEW_HYPS(i).COST = roomAlignmentMex(HYP_COMP, GND_COMP, 3, 0);
            end 
        end
        
    end
    
end

save([bufname config.localResultFile], 'ITERHYPS');


end

