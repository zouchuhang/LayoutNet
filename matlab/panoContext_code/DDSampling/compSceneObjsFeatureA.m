function [ sceneImgFea ] = compSceneObjsFeatureA( ALL_HYPS, OBJFEA )
%COMPSCENEHYPSFEATURE Summary of this function goes here
%   Detailed explanation goes here
num_hyps = length(ALL_HYPS);
% coor = [];
% load('./rectangleDetector/segmentation/uniformvector_lvl6.mat');
% [omH, omW, ~] = size(omap);
% coor2D = uv2coords(xyz2uvN(coor), omW, omH);
% coor2D_indi = sub2ind([omH omW], coor2D(:,2), coor2D(:,1));
% vector_num = size(coor2D_indi,1);

sceneImgFea = zeros(num_hyps, 2830);

% OBJFEA = hyps_obj.objfea;
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

roomdescriptor = zeros(1,5);
numberfeature = zeros(1,37);

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

    for i = 1:num_type
        if ~sel_hyps(i).fixed
            continue;
        end
        v = sel_hyps(i).selID;
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
            [sum(gc_sum_local,1) mean(gc_sum_local, 1) max(gc_sum_local,[],1) min(gc_sum_local,[],1) ...
             sum(gc_avg_local,1) mean(gc_avg_local, 1) max(gc_avg_local,[],1) min(gc_avg_local,[],1) ...
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
        cn_feature(i,:) = [sum(cn_entropy_local) mean(cn_entropy_local) max(cn_entropy_local) min(cn_entropy_local)];
        % for every object, we take its area
        indi_local = sz(v);
        sz_feature(i,:) = [sum(indi_local) mean(indi_local) max(indi_local) min(indi_local)];        
    end
    
        
    gc_global_feature = zeros(2,33);
    om_global_feature = zeros(2,15);
    cn_global_feature = zeros(2,1);
    sz_global_feature = zeros(2,1);
    
%     gc_global_feature(1,:) = [gc_global_fore_sum(hid,:) gc_global_fore_avg(hid,:) gc_global_fore_std(hid,:)];
%     gc_global_feature(2,:) = [gc_global_back_sum(hid,:) gc_global_back_avg(hid,:) gc_global_back_std(hid,:)];
%     om_global_feature(1,:) = [om_global_fore_sum(hid,:) om_global_fore_avg(hid,:) om_global_fore_std(hid,:)];
%     om_global_feature(2,:) = [om_global_back_sum(hid,:) om_global_back_avg(hid,:) om_global_back_std(hid,:)];
%     cn_global_feature = [cn_entropy_fore_local(hid);cn_entropy_back_local(hid)];
%     sz_global_feature = [sz_fore_global(hid);sz_back_global(hid)];
    
    
    
    sceneImgFea(hid,:) = [reshape([sm_feature bu_feature gc_feature om_feature cn_feature sz_feature]',1,[]) ...
                          reshape([gc_global_feature om_global_feature cn_global_feature sz_global_feature]',1,[]) ...
                          roomdescriptor numberfeature];
       
end

end
%%
function p = sigmoidFunc_flex(v, a, b)
% p = 1./(1+exp(a*v+b));
p = 1./(1+exp(a*(v-b)));
end
