folderName = 'annotation/dataset/bedroom/'; % buffer dataset
load([folderName 'ANNO_ALL.mat']);
load([folderName 'IMGLIST.mat']);
data_valid = true(length(ANNO_ALL), 1);
for i = 1:length(ANNO_ALL)
%     fprintf('%d\n',i);
    if isempty(ANNO_ALL(i).ANNO3D) || ~ANNO_ALL(i).ANNO3D.b_singleroom
        data_valid(i) = false;
    else
        data_valid(i) = true;
    end
end
VALID_SAMPLE_ID = find(data_valid);
valid_sample_number = length(VALID_SAMPLE_ID);

%%
try load('fullfea_t.mat')
catch
% total_fea_dim = comp_object_3d_feature() + comp_object_img_feature();
total_fea_dim = comp_object_3d_feature();
DATA = repmat(struct('X',zeros(0,total_fea_dim),'Y',cell(0,1),'T',zeros(0,1),'aid',[],'angle',[]), valid_sample_number, 1);
for j = 1:valid_sample_number
    fprintf('%d/%d\n', j, valid_sample_number);
    aid = VALID_SAMPLE_ID(j);
    DATA(j).aid = aid;
    if isempty(ANNO_ALL(aid).ANNO3D) || ~ANNO_ALL(aid).ANNO3D.b_singleroom
        continue;
    end
    objects3D = ANNO_ALL(aid).ANNO3D.objects3D;
    contain3Dcoords = ANNO_ALL(aid).ANNO3D.contain3Dcoords;
    room_id = contain3Dcoords(1);
    room_points = objects3D(room_id).out_points_w;
    room_xyz = [min(room_points, [], 1) max(room_points, [],  1)];
    
%     % decide any rotation
%     bed_id = [];
%     for i = 2:length(contain3Dcoords)
%         if strcmp('bed',objects3D(contain3Dcoords(i)).name)
%             bed_id = contain3Dcoords(i);
%             break;
%         end
%     end
%     if ~isempty(bed_id)
%         bed_points = objects3D(bed_id).out_points_w;
%         bed_xyz = [min(bed_points, [], 1) max(bed_points, [],  1)];
%         dist_diff = bed_xyz - room_xyz;
%         [~, id] = min(abs(dist_diff([1 2 4 5])));
%         angle = [0 pi pi/2 -pi/2];
%         angle = angle(id);
%         for i = 1:length(contain3Dcoords)
%             p = objects3D(contain3Dcoords(i)).out_points_w;
%             r = p;
%             r(:,1) = p(:,1).*cos(angle) + p(:,2)*sin(angle);
%             r(:,2) = -p(:,1).*sin(angle) + p(:,2)*cos(angle);
%             objects3D(contain3Dcoords(i)).out_points_w = r;
%         end
%         DATA(j).angle = angle;
%     else
%         DATA(j).angle = 0;
%     end       
%     room_points = objects3D(room_id).out_points_w;
%     room_xyz = [min(room_points, [], 1) max(room_points, [],  1)];
%     % end
    
    DATA(j).X = zeros(length(contain3Dcoords)-1,total_fea_dim);
    DATA(j).Y = cell(length(contain3Dcoords)-1,1);
    DATA(j).T = zeros(length(contain3Dcoords)-1,1);
    for i = 2:length(contain3Dcoords)
        oid = contain3Dcoords(i);
%         if (objects3D(oid).type==0 || objects3D(oid).type==1)
%             continue;
%         end
        DATA(j).T(i-1,1) = objects3D(oid).type;
        obj_points = objects3D(oid).out_points_w;
        obj_xyz = [min(obj_points, [], 1) max(obj_points, [], 1)];
        obj_x = objects3D(oid).x_w;
        obj_ori_size = obj_x((i-1)*7+4:(i-1)*7+6)';
        obj_angle = obj_x((i-1)*7+7);
        feature = comp_object_3d_feature(obj_xyz, room_xyz, obj_ori_size, obj_angle);              
%         img_feature = comp_object_img_feature(img, ANNO_ALL(aid).anno3D(oid).points, ANNO_ALL(aid).anno3D(oid).type);
        
%         DATA(j).X(i-1,:) = [feature img_feature];
        DATA(j).X(i-1,:) = feature;
        DATA(j).Y{i-1,1} = objects3D(oid).name;
    end
end
end

%% statistics
all_hist = zeros(28,100);
for i = 1:valid_sample_number
    label = get_object_type( DATA(i).Y );
    h = hist(label, 1:28);
    for j = 1:28
        all_hist(j,h(j)+1) = all_hist(j,h(j)+1)+1;
    end

end


%%
TRAIN_SAMPLE_IND = randsample(valid_sample_number, round(valid_sample_number*0.8))';
TEST_SAMPLE_IND = setdiff(1:valid_sample_number, TRAIN_SAMPLE_IND);

TRAINX = zeros(0,total_fea_dim);
TRAINY = cell(0,1);
for i = TRAIN_SAMPLE_IND
    fprintf('%d\n',i);
%     valid = ~(DATA(i).T==0 | DATA(i).T==1);
    valid = true(length(DATA(i).Y),1);
%     num = length(DATA(i).Y);
    num = sum(valid);
    TRAINX(end+1:end+num,:) = DATA(i).X(valid,:);
    TRAINY(end+1:end+num,1) = DATA(i).Y(valid);
end
TRAIN_LABEL = get_object_type( TRAINY );
% TRAIN_LABEL(TRAIN_LABEL>=13) = 28;
TESTX = zeros(0,total_fea_dim);
TESTY = cell(0,1);
for i = TEST_SAMPLE_IND
    fprintf('%d\n',i);
%     valid = ~(DATA(i).T==0 | DATA(i).T==1);
    valid = true(length(DATA(i).Y),1);
%     num = length(DATA(i).Y);
    num = sum(valid);
    TESTX(end+1:end+num,:) = DATA(i).X(valid,:);
    TESTY(end+1:end+num,1) = DATA(i).Y(valid);
end
TEST_LABEL = get_object_type( TESTY );
% TEST_LABEL(TEST_LABEL>=13) = 28;

%% training
% B = TreeBagger(100,TRAINX,TRAIN_LABEL,'Method','classification','FBoot',0.6,'NVarToSample',7);
% B = TreeBagger(100,TRAINX,TRAIN_LABEL,'Method','classification','FBoot',0.6);
% train_err = error(B, TRAINX, TRAIN_LABEL, 'mode', 'cumulative');
% test_err = error(B, TESTX, TEST_LABEL, 'mode', 'cumulative');
% B = TreeBagger(100,[TRAINX;TESTX],[TRAIN_LABEL;TEST_LABEL],'Method','classification','FBoot',0.6, 'OOBPred', 'on');
B = TreeBagger(100,[TRAINX],[TRAIN_LABEL],'Method','classification','FBoot',0.6, 'OOBPred', 'on');
err = oobError(B);
% TRAINX = TRAINX(:,1:65);
% TESTX = TESTX(:,1:65);

% B_all = TreeBagger(100,TRAINX,TRAIN_LABEL,'Method','classification','FBoot',0.6);
% TRAIN_LABEL_MERGE = TRAIN_LABEL;
% TRAIN_LABEL_MERGE(TRAIN_LABEL_MERGE>=13) = 28;
% B_mgr = TreeBagger(100,TRAINX,TRAIN_LABEL_MERGE,'Method','classification','FBoot',0.6);
% [ test_err_all, ambimat_all, LABELERROR_all, LABELNUM_all, fea_vote_all, d_all ] = result_analysis( TESTX, TEST_LABEL, B_all );
% TEST_LABEL_MERGE = TEST_LABEL;
% TEST_LABEL_MERGE(TEST_LABEL_MERGE>=13) = 28;
% [ test_err_mgr, ambimat_mgr, LABELERROR_mgr, LABELNUM_mgr, fea_vote_mgr, d_mgr ] = result_analysis( TESTX, TEST_LABEL_MERGE, B_mgr );
% figure(1); clf; 
% bar(1:classnum, [LABELNUM_all;LABELNUM_all-LABELERROR_all;LABELNUM_mgr-LABELERROR_mgr]'); 
% title(sprintf('All TYPE: %f, MERGE TYPE: %f', test_err_all(end), test_err_mgr(end)));
% set(gca, 'XTick', 1:classnum);
% LABELSTRING = get_object_type([1:classnum]);
% set(gca, 'XTickLabel', LABELSTRING);
% xticklabel_rotate([],90,[],'Fontsize',11);


% B_full = TreeBagger(100,TRAINX,TRAIN_LABEL,'Method','classification','FBoot',0.6, 'OOBPred', 'on');
% err_full = oobError(B_full);
% B_sktk = TreeBagger(100,TRAINX(:,1:65),TRAIN_LABEL,'Method','classification','FBoot',0.6, 'OOBPred', 'on');
% err_sktk = oobError(B_sktk);
% B_3dol = TreeBagger(100,TRAINX(:,1:55),TRAIN_LABEL,'Method','classification','FBoot',0.6, 'OOBPred', 'on');
% err_3dol = oobError(B_3dol);
% 
% [ test_err_full, ambimat_full, LABELERROR_full, LABELNUM_full, fea_vote_full, d_full ] = result_analysis( TESTX, TEST_LABEL, B_full );
% [ test_err_sktk, ambimat_sktk, LABELERROR_sktk, LABELNUM_sktk, fea_vote_sktk, d_sktk ] = result_analysis( TESTX(:,1:65), TEST_LABEL, B_sktk );
% [ test_err_3dol, ambimat_3dol, LABELERROR_3dol, LABELNUM_3dol, fea_vote_3dol, d_3dol ] = result_analysis( TESTX(:,1:55), TEST_LABEL, B_3dol );
% 
% figure(1); clf; 
% bar(1:classnum, [LABELNUM_full;LABELNUM_full-LABELERROR_full;LABELNUM_sktk-LABELERROR_sktk;LABELNUM_3dol-LABELERROR_3dol]'); 
% title(sprintf('3D Error: %f, +SketchToken: %f, +SIFT: %f', test_err_3dol(end), test_err_sktk(end), test_err_full(end)));
% set(gca, 'XTick', 1:classnum);
% LABELSTRING = get_object_type([1:classnum]);
% set(gca, 'XTickLabel', LABELSTRING);
% xticklabel_rotate([],90,[],'Fontsize',11);



%% 
% confusion matrix
testgt = TEST_LABEL;
test_err = error(B, TESTX, TEST_LABEL, 'mode', 'cumulative');
[PREDICTY,scores,stdevs] = predict(B,TESTX);
[~,testre] = max(scores, [], 2);
ClassNames = B.ClassNames;
for i = 1:length(testre)
    testre(i) = str2double(ClassNames{testre(i)});
end

classnum = 12;
ambimat = zeros(classnum, classnum);
for i = 1:length(TESTY)
    if testgt(i)>classnum || testre(i)>classnum
        continue;
    end
    if testgt(i) ~= testre(i);
        ambimat(testgt(i), testre(i)) = ambimat(testgt(i), testre(i)) + 1;
    end
end
%%
type_num = hist(testgt, 1:classnum);
norm_ambimat = ambimat./repmat(type_num'+0.000000001, 1, classnum);
figure(1); clf; imagesc(norm_ambimat);
set(gca, 'XTick', 1:classnum);
set(gca, 'YTick', 1:classnum);
LABELSTRING = get_object_type([1:classnum]);
set(gca, 'XTickLabel', LABELSTRING);
set(gca, 'YTickLabel', LABELSTRING);
xticklabel_rotate([],90,[],'Fontsize',11);

%% accuracy
LABELERROR = sum(ambimat,2)';
LABELNUM = hist(testgt, 1:classnum);
figure(2); clf; bar(1:classnum, [LABELNUM;LABELNUM-LABELERROR]'); title(sprintf('Error: %f', test_err(end)));
set(gca, 'XTick', 1:classnum);
set(gca, 'XTickLabel', LABELSTRING);
xticklabel_rotate([],90,[],'Fontsize',11);

[~,I] = max(norm_ambimat,[],2);
templabel = get_object_type(I);
d = [LABELSTRING templabel];

% % feature vote
% feadim = size(TRAINX,2);
% Trees = B.Trees;
% fea_vote = zeros(feadim,1);
% for tid = 1:length(Trees)
%     T = Trees{tid};
%     v = cutvar(T);
%     for lid = 1:200%length(v)
%         if ~isempty(v{lid})
%             n = str2double(v{lid}(2:end));
%             fea_vote(n) = fea_vote(n) + 1;
%         end
%     end
% end
% figure(3); clf; bar(1:feadim, fea_vote);