function [objfea, roomfea, unaryanglesid] = compObjHypsFeature( rotImg, omap, gc, cn, CAN_POOL, typenum, bufName)
%COMPOBJHYPSFEATURE Summary of this function goes here
%   Detailed explanation goes here
% if ~exist('unaryname','var')
%     unaryname = './object_classifier/objectunary.mat';
% end
global config;

all_hyps = CAN_POOL.sel_hyps;
room_obj = CAN_POOL.room;

num_object = length(all_hyps.objects);
objfea = repmat(struct('rectscore',[],'segmscore',[], ...
    'randomforest',[],'unaryscore',[], ...
    'indicator',[],'om_sum',[],'om_avg', [], 'om_std', [], ...
    'gc_sum',[],'gc_avg', [], 'gc_std', [], ...
    'cn_hist', [], 'cn_entropy', [], 'sz',[]), ...
    num_object, 1);
room_xyz = [min(room_obj.out_points_w,[],1) max(room_obj.out_points_w,[],1)];

%% get image score
checkrect = true(1,1);
checkseg = true(1,1);
objects = all_hyps.objects;
[ cuboid_rect_score, segment_rect_score] = ...
    computeCuboidImgFea( rotImg, {[objects(1); objects]}, checkrect, checkseg, bufName );

for oid = 1:num_object
    objfea(oid).rectscore = max(cuboid_rect_score{1}(oid+1).surfscore);
    if isempty(objfea(oid).rectscore)
        objfea(oid).rectscore = -10;
    end
    objfea(oid).segmscore = segment_rect_score{1}(oid+1);
end

%% get semantic score
% load('./object_classifier/objectunary.mat');
load(config.unaryName);
[ unaryscore, unaryanglesid ] = unarySizeScore( room_xyz, all_hyps.obj_xyz, typenum, unary_size);
forcecuboid = [2 7 9 10 12];
forcerect = [1 3 4 8];
isrectangle = any((all_hyps.obj_xyz(:,4:6)-all_hyps.obj_xyz(:,1:3))<0.01,2);
iscuboid = ~isrectangle;
unaryscore(isrectangle, forcecuboid) = 1000000;
unaryscore(iscuboid, forcerect) = 1000000;

for oid = 1:num_object
    objfea(oid).unaryscore = unaryscore(oid,:);
    objfea(oid).randomforest = all_hyps.scores(oid,:);
end

roomfea.unaryscore = unarySizeScore( room_xyz, [], 0, unary_size);
roomfea.omapScr = room_obj.omapScr;
roomfea.gcScr = room_obj.gcScr;
roomfea.mgScr = room_obj.mgScr;
roomfea.opScr = room_obj.opScr;

%% get vectorized marker
coor = [];
% load('./rectangleDetector/segmentation/uniformvector_lvl6.mat');
[ coor, tri ] = getUniformVector( 6 );
[ indicator ] = insideIndicator( all_hyps.objects, coor );
for oid = 1:num_object
    objfea(oid).indicator = indicator(:,oid);
end

%% get omap and gc feature
[omH, omW, C] = size(omap);
if C~=3
    fprintf('Warning: OMAP should have 3 channels.\n');
end
coor2D = uv2coords(xyz2uvN(coor), omW, omH);
coor2D_indi = sub2ind([omH omW], coor2D(:,2), coor2D(:,1));
[ om_sum, om_avg, om_std ] = omapFeature( omap, indicator, coor2D_indi );
[ gc_sum, gc_avg, gc_std  ] = gcFeature( gc, indicator, coor2D_indi );
[ cn_hist, cn_entropy ] = colorNameFeature( cn, indicator, coor2D_indi );
sz_vec = sum(indicator,1)';
for oid = 1:num_object
    objfea(oid).om_sum = om_sum(oid,:);
    objfea(oid).om_avg = om_avg(oid,:);
    objfea(oid).om_std = om_std(oid,:);
    
    objfea(oid).gc_sum = gc_sum(oid,:);
    objfea(oid).gc_avg = gc_avg(oid,:);
    objfea(oid).gc_std = gc_std(oid,:);
    
    objfea(oid).cn_hist = cn_hist(oid,:);
    objfea(oid).cn_entropy = cn_entropy(oid);
    objfea(oid).sz = sz_vec(oid);
end

%% 

end

