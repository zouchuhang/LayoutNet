function [ ALL_HYPS, CAN_POOL ] = groundTruthHypothesisNewFunc( gnd, omap, gc, cn, rotImg, typenum, bufName)
%GROUNDTRUTHHYPOTHESISNEWFUNC Summary of this function goes here
%   Detailed explanation goes here
global config;
%% compute the checking direction on omap                     
% [candiSetXYZ, ~] = icosahedron2sphere(6);
[candiSetXYZ, ~] = getUniformVector( 6 );
[ohei, owid, ~] = size(omap); 
candiSetUV = uv2coords(xyz2uvN(candiSetXYZ), owid, ohei);
candiInd = sub2ind([ohei owid], candiSetUV(:,2), candiSetUV(:,1));
% convert gc to omap
gndGC = zeros(ohei, owid, 3);
gndGC(:,:,1) = max(gc(:,:,5),gc(:,:,6));
gndGC(:,:,2) = max(gc(:,:,1),gc(:,:,3));
gndGC(:,:,3) = max(gc(:,:,2),gc(:,:,4));
normGndGC = sum(gndGC(candiInd))+sum(gndGC(candiInd+1*ohei*owid))+sum(gndGC(candiInd+2*ohei*owid));
% omap
gndOmap = omap;
gndOmap(gndOmap>0) = 1;
normGndOmap = sum(gndOmap, 3);
gndOmap = gndOmap./(repmat(normGndOmap+0.0001, [1 1 3])); 
normGndOmap = sum(gndOmap(candiInd))+sum(gndOmap(candiInd+1*ohei*owid))+sum(gndOmap(candiInd+2*ohei*owid));
numPoint = size(candiSetXYZ,1);

% merge
gndMerge = zeros(ohei, owid, 3);
gndMerge(1:ohei/2,:,:) = gndOmap(1:ohei/2,:,:);
gndMerge(ohei/2+1:ohei, :, :) = gndGC(ohei/2+1:ohei, :, :);
normMegMap = sum(gndMerge(candiInd))+sum(gndMerge(candiInd+1*ohei*owid))+sum(gndMerge(candiInd+2*ohei*owid));

% optimal merge
gndOptim = zeros(ohei, owid, 3);
gndOptim(1:685,:,:) = gndOmap(1:685,:,:);
gndOptim(686:ohei, :, :) = gndGC(686:ohei, :, :);
normOptMap = sum(gndOptim(candiInd))+sum(gndOptim(candiInd+1*ohei*owid))+sum(gndOptim(candiInd+2*ohei*owid));

%%
room = gnd(1);
% compute room hyp omap feature
hypVp = false(numPoint,3);
align = room.align;
valid1 = insideCone(align([6;5;8;7],:), candiSetXYZ, 0);
valid2 = insideCone(align([7;8;4;3],:), candiSetXYZ, 0);
valid3 = insideCone(align([6;7;3;2],:), candiSetXYZ, 0);
valid4 = insideCone(align([3;4;1;2],:), candiSetXYZ, 0);
valid5 = insideCone(align([5;6;2;1],:), candiSetXYZ, 0);
valid6 = insideCone(align([8;5;1;4],:), candiSetXYZ, 0);
hypVp(valid1,1) = true; hypVp(valid4,1) = true;
hypVp(valid2,2) = true; hypVp(valid5,2) = true;
hypVp(valid3,3) = true; hypVp(valid6,3) = true;
response = sum(gndOmap(candiInd(hypVp(:,1))+0*ohei*owid)) ...
         + sum(gndOmap(candiInd(hypVp(:,2))+1*ohei*owid)) ...
         + sum(gndOmap(candiInd(hypVp(:,3))+2*ohei*owid));
room.omapScr = response/normGndOmap;
response = sum(gndGC(candiInd(hypVp(:,1))+0*ohei*owid)) ...
         + sum(gndGC(candiInd(hypVp(:,2))+1*ohei*owid)) ...
         + sum(gndGC(candiInd(hypVp(:,3))+2*ohei*owid));
room.gcScr = response/normGndGC;
response = sum(gndMerge(candiInd(hypVp(:,1))+0*ohei*owid)) ...
         + sum(gndMerge(candiInd(hypVp(:,2))+1*ohei*owid)) ...
         + sum(gndMerge(candiInd(hypVp(:,3))+2*ohei*owid));
room.mgScr = response/normMegMap;
response = sum(gndOptim(candiInd(hypVp(:,1))+0*ohei*owid)) ...
         + sum(gndOptim(candiInd(hypVp(:,2))+1*ohei*owid)) ...
         + sum(gndOptim(candiInd(hypVp(:,3))+2*ohei*owid));
room.opScr = response/normOptMap;

%% 
% load('./object_classifier/objectclassifier.mat');
load(config.classifierName);
objects = gnd(2:end);
[ scores, valid ] = objectClassifier(room, objects, B);
obj_xyz = zeros(length(objects),6);
for i = 1:length(objects)
    obj_xyz(i,:) = [objects(i).align(1,:) objects(i).align(7,:)];
end
sel_hyps.objects = objects;
sel_hyps.scores = scores;
sel_hyps.valid = valid;
sel_hyps.obj_xyz = obj_xyz;

CAN_POOL{1}.sel_hyps = sel_hyps;
CAN_POOL{1}.room = room;

sel_hyps = repmat(struct('fixed',false,'selID',[]), typenum, 1);
for i = 1:length(objects)
    if objects(i).objtype>12
        continue;
    end
    sel_hyps(objects(i).objtype).fixed = true;
    sel_hyps(objects(i).objtype).selID = [sel_hyps(objects(i).objtype).selID i];
end
ALL_HYPS(1).sel_hyps = sel_hyps;


[objfea, roomfea, anglesid] = compObjHypsFeature( rotImg, omap, gc, cn, CAN_POOL{1}, typenum, bufName );
CAN_POOL{1}.sel_hyps.objfea = objfea;
CAN_POOL{1}.room.roomfea = roomfea;
CAN_POOL{1}.sel_hyps.anglesid = anglesid;

[ sceneImgFea ] = compSceneHypsFeatureA( ALL_HYPS, CAN_POOL{1}.sel_hyps, CAN_POOL{1}.room, omap, gc, cn );
ALL_HYPS(1).sceneImgFea = sceneImgFea;

end

