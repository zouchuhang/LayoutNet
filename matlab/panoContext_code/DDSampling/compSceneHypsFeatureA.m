function [ sceneImgFea ] = compSceneHypsFeatureA( ALL_HYPS, hyps_obj, room_obj, omap, gc, cn )
%COMPSCENEHYPSFEATURE Summary of this function goes here
%   Detailed explanation goes here
global config

num_hyps = length(ALL_HYPS);
% coor = [];
% load('./rectangleDetector/segmentation/uniformvector_lvl6.mat');
[ coor, tri ] = getUniformVector( 6 );
[omH, omW, ~] = size(omap);
coor2D = uv2coords(xyz2uvN(coor), omW, omH);
coor2D_indi = sub2ind([omH omW], coor2D(:,2), coor2D(:,1));
vector_num = size(coor2D_indi,1);

% sceneImgFea = zeros(num_hyps, 2271);
sceneImgFea = zeros(num_hyps, 2830);

OBJFEA = hyps_obj.objfea;
unaryscore = vertcat(OBJFEA.unaryscore);
rectscore = vertcat(OBJFEA.rectscore);
segmscore = vertcat(OBJFEA.segmscore);
randomforest = vertcat(OBJFEA.randomforest);

unaryscore = sigmoidFunc_flex(unaryscore, 33.67, 0.15);
rectscore = sigmoidFunc_flex(rectscore, -4, 0);
% segmscore = ones(length(segmscore),1);
segmscore = sigmoidFunc_flex(segmscore, -4, 0.5);
randomforest = sigmoidFunc_flex(randomforest, -4, 0.5);


om_sum = vertcat(OBJFEA.om_sum);
om_avg = vertcat(OBJFEA.om_avg);
om_std = vertcat(OBJFEA.om_std);
gc_sum = vertcat(OBJFEA.gc_sum);
gc_avg = vertcat(OBJFEA.gc_avg);
gc_std = vertcat(OBJFEA.gc_std);
cn_entropy = vertcat(OBJFEA.cn_entropy);
sz = vertcat(OBJFEA.sz);

% bu = rectscore.*segmscore;
% sm = unaryscore(:,1:12).*randomforest(:,1:12);
indicator = horzcat(OBJFEA.indicator);

%% compute global feature
indicator_global = false(vector_num, num_hyps);
for hid = 1:num_hyps
    sel_hyps = ALL_HYPS(hid).sel_hyps;
    v = horzcat(sel_hyps.selID);
    indicator_global(:,hid) = indicator_global(:,hid) | any(indicator(:,v),2);
%     num_type = length(sel_hyps);
%     for i = 1:num_type
%         if ~sel_hyps(i).fixed
%             continue;
%         end
%         v = sel_hyps(i).selID;
% %         fprintf('%d ',v);
%         indicator_global(:,hid) = indicator_global(:,hid) | any(indicator(:,v),2);
%     end
end
[ om_global_fore_sum, om_global_fore_avg, om_global_fore_std ] = omapFeature( omap, indicator_global, coor2D_indi );
[ gc_global_fore_sum, gc_global_fore_avg, gc_global_fore_std  ] = gcFeature( gc, indicator_global, coor2D_indi );
[ ~, cn_entropy_fore_local ] = colorNameFeature( cn, indicator_global, coor2D_indi );
sz_fore_global = sum(indicator_global,1)';

[ om_global_back_sum, om_global_back_avg, om_global_back_std ] = omapFeature( omap, ~indicator_global, coor2D_indi );
[ gc_global_back_sum, gc_global_back_avg, gc_global_back_std  ] = gcFeature( gc, ~indicator_global, coor2D_indi );
[ ~, cn_entropy_back_local ] = colorNameFeature( cn, ~indicator_global, coor2D_indi );
sz_back_global = sum(~indicator_global,1)';

% indicator_wall_room = any(room_obj.roomfea.indicator(:,[2 3 5 6]),2);
% indicator_floor_room = room_obj.roomfea.indicator(:,4);
% indicator_ceil_room = room_obj.roomfea.indicator(:,1);
% indicator_wall = repmat(indicator_wall_room, 1, num_hyps) & ~indicator_global;
% indicator_floor = repmat(indicator_floor_room, 1, num_hyps) & ~indicator_global;
% indicator_ceil = repmat(indicator_ceil_room, 1, num_hyps) & ~indicator_global;
% 
% [ om_global_wall_sum, om_global_wall_avg, om_global_wall_std ] = omapFeature( omap, indicator_wall, coor2D_indi );
% [ gc_global_wall_sum, gc_global_wall_avg, gc_global_wall_std  ] = gcFeature( gc, indicator_wall, coor2D_indi );
% [ ~, cn_entropy_wall_local ] = colorNameFeature( cn, indicator_wall, coor2D_indi );
% sz_wall_global = sum(indicator_wall,1)';
% 
% [ om_global_floor_sum, om_global_floor_avg, om_global_floor_std ] = omapFeature( omap, indicator_floor, coor2D_indi );
% [ gc_global_floor_sum, gc_global_floor_avg, gc_global_floor_std  ] = gcFeature( gc, indicator_floor, coor2D_indi );
% [ ~, cn_entropy_floor_local ] = colorNameFeature( cn, indicator_floor, coor2D_indi );
% sz_floor_global = sum(indicator_floor,1)';
% 
% [ om_global_ceil_sum, om_global_ceil_avg, om_global_ceil_std ] = omapFeature( omap, indicator_ceil, coor2D_indi );
% [ gc_global_ceil_sum, gc_global_ceil_avg, gc_global_ceil_std  ] = gcFeature( gc, indicator_ceil, coor2D_indi );
% [ ~, cn_entropy_ceil_local ] = colorNameFeature( cn, indicator_ceil, coor2D_indi );
% sz_ceil_global = sum(indicator_ceil,1)';

%% short room descriptor

% roomfea.unaryscore = unarySizeScore( room_xyz, [], 0, unary_size);
% roomfea.omapScr = room_obj.omapScr;
% roomfea.gcScr = room_obj.gcScr;
% roomfea.mgScr = room_obj.mgScr;
% roomfea.opScr = room_obj.opScr;

ROOMFEA = room_obj.roomfea;
roomdescriptor = [ROOMFEA.unaryscore ROOMFEA.omapScr ROOMFEA.gcScr ROOMFEA.mgScr ROOMFEA.opScr];
% load('./object_classifier/objectnumber.mat');
load(config.numberName);
maxprob = max(prob,[],2);


%% objects
for hid = 1:num_hyps
    sel_hyps = ALL_HYPS(hid).sel_hyps;
    num_type = length(sel_hyps);
    
    sm_feature = zeros(num_type,12);
    bu_feature = zeros(num_type,12);
    gc_feature = zeros(num_type,132);
    om_feature = zeros(num_type,60);
    cn_feature = zeros(num_type,4);
    sz_feature = zeros(num_type,4);
    nb_feature = zeros(num_type,1);
    
    
    for i = 1:num_type
        if ~sel_hyps(i).fixed
            continue;
        end
        v = sel_hyps(i).selID;
        nb_feature(i) = length(v);
        sm = [unaryscore(v,i) randomforest(v,i)];
        bu = [rectscore(v) segmscore(v)];
        
        sm_local = [sm prod(sm,2)];
        sm_feature(i,:) = [sum(sm_local,1) mean(sm_local,1) max(sm_local,[],1) min(sm_local,[],1)];
        bu_local = [bu prod(bu,2)];
        bu_feature(i,:) = [sum(bu_local,1) mean(bu_local,1) max(bu_local,[],1) min(bu_local,[],1)];
        % for every object, we take gc and sum(1:4) sum(1:6) as feature 
        gc_sum_local = gc_sum(v,:);
        gc_avg_local = gc_avg(v,:);
        gc_std_local = gc_std(v,:);
        gc_feature(i,:) = ...
            [sum(gc_sum_local,1) mean(gc_sum_local, 1) max(gc_sum_local,[],1) min(gc_sum_local,[],1)...
             sum(gc_avg_local,1) mean(gc_avg_local, 1) max(gc_avg_local,[],1) min(gc_avg_local,[],1)...
             sum(gc_std_local,1) mean(gc_std_local, 1) max(gc_std_local,[],1) min(gc_std_local,[],1)]; 
        % for every object, we take om and sum(1:2)
        om_sum_local = om_sum(v,:);
        om_avg_local = om_avg(v,:);
        om_std_local = om_std(v,:);
        om_feature(i,:) = ...
            [sum(om_sum_local,1) mean(om_sum_local, 1) max(om_sum_local,[],1) min(om_sum_local,[],1)...
             sum(om_avg_local,1) mean(om_avg_local, 1) max(om_avg_local,[],1) min(om_avg_local,[],1)...
             sum(om_std_local,1) mean(om_std_local, 1) max(om_std_local,[],1) min(om_std_local,[],1)]; 
        % for every object, we take entropy of color name histogram
        cn_entropy_local = cn_entropy(v); %-sum(cn_hist(v,:).*log2(cn_hist(v,:)+0.0001), 2);
        cn_feature(i,:) = [sum(cn_entropy_local) mean(cn_entropy_local) ...
            max(cn_entropy_local) min(cn_entropy_local)];
        % for every object, we take its area
        indi_local = sz(v);
        sz_feature(i,:) = [sum(indi_local) mean(indi_local) max(indi_local) min(indi_local)];        
    end
    
        
    gc_global_feature = zeros(2,33);
    om_global_feature = zeros(2,15);
%     cn_global_feature = zeros(2,1);
%     sz_global_feature = zeros(2,1);
    
    gc_global_feature(1,:) = [gc_global_fore_sum(hid,:) gc_global_fore_avg(hid,:) gc_global_fore_std(hid,:)];
    gc_global_feature(2,:) = [gc_global_back_sum(hid,:) gc_global_back_avg(hid,:) gc_global_back_std(hid,:)];
    om_global_feature(1,:) = [om_global_fore_sum(hid,:) om_global_fore_avg(hid,:) om_global_fore_std(hid,:)];
    om_global_feature(2,:) = [om_global_back_sum(hid,:) om_global_back_avg(hid,:) om_global_back_std(hid,:)];
    cn_global_feature = [cn_entropy_fore_local(hid);cn_entropy_back_local(hid)];
    sz_global_feature = [sz_fore_global(hid);sz_back_global(hid)];
    
    nb_prior_feature = zeros(1, num_type);
    nb_rela_prior_feature = zeros(1, num_type);
    for i = 1:num_type
        nb_prior_feature(i) = prob(i,nb_feature(i)+1);
        nb_rela_prior_feature(i) = prob(i,nb_feature(i)+1) / maxprob(i);
    end
    
    
    
    
    sceneImgFea(hid,:) = [reshape([sm_feature bu_feature gc_feature om_feature cn_feature sz_feature]',1,[]) ...
                      reshape([gc_global_feature om_global_feature cn_global_feature sz_global_feature]',1,[]) ...
                      roomdescriptor nb_feature' sum(nb_feature) nb_prior_feature nb_rela_prior_feature];
    
%     gc_room_feature = zeros(3,33);
%     om_room_feature = zeros(3,15);
%     gc_room_feature(1,:) = [gc_global_wall_sum(hid,:) gc_global_wall_avg(hid,:) gc_global_wall_std(hid,:)];
%     gc_room_feature(2,:) = [gc_global_floor_sum(hid,:) gc_global_floor_avg(hid,:) gc_global_floor_std(hid,:)];
%     gc_room_feature(3,:) = [gc_global_ceil_sum(hid,:) gc_global_ceil_avg(hid,:) gc_global_ceil_std(hid,:)];
%     om_room_feature(1,:) = [om_global_wall_sum(hid,:) om_global_wall_avg(hid,:) om_global_wall_std(hid,:)];
%     om_room_feature(2,:) = [om_global_floor_sum(hid,:) om_global_floor_avg(hid,:) om_global_floor_std(hid,:)];
%     om_room_feature(3,:) = [om_global_ceil_sum(hid,:) om_global_ceil_avg(hid,:) om_global_ceil_std(hid,:)];
%     cn_room_feature = [cn_entropy_wall_local(hid); cn_entropy_floor_local(hid); cn_entropy_ceil_local];
%     sz_room_feature = [sz_wall_global(hid); sz_floor_global(hid); sz_ceil_global(hid)];
        
%     sceneImgFea(hid,:) = [reshape([sm_feature bu_feature gc_feature om_feature cn_feature sz_feature]',1,[]) ...
%                           reshape([gc_global_feature om_global_feature cn_global_feature sz_global_feature]',1,[]) ...
%                           reshape([gc_room_feature om_room_feature cn_room_feature sz_room_feature]',1,[]) ...
%                           roomdescriptor];

       
end

end
%%
function p = sigmoidFunc_flex(v, a, b)
% p = 1./(1+exp(a*v+b));
p = 1./(1+exp(a*(v-b)));
end
