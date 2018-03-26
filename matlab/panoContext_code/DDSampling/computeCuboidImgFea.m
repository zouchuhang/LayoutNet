function [ cuboid_rect_score, segment_rect_score] = computeCuboidImgFea( rotImg, allhyps, checkrect, checkseg, bufName )
%COMPUTECUBOIDIMGFEA Summary of this function goes here
%   Detailed explanation goes here

%% project image to 6 views and compute pixel-wise hog
rotImg  = im2double(rotImg);
vp = [-1  0  0; ...
       1  0  0; ...
       0 -1  0; ...
       0  1  0; ...
       0  0 -1; ...
       0  0  1];
uv = xyz2uvN(vp);
fov = 160/180*pi;
cutSize = 2000;

hasBuff = exist('bufName','var') && exist(bufName,'file') && ~isempty(bufName);
% hasBuff = false;
if hasBuff
    load(bufName);
else
    all_hog = cell(6,1);
    [sepSceneSix] = separatePano( rotImg, fov, uv(:,1), uv(:,2), cutSize);
end
% load './rectangleDetector/finalModel.mat'; % load in allModel
% rectModel = allModel{4}; % after 2nd negative mining
% load './rectangleDetector/selectModel.mat';

% compute response
detconfig;
load(config.modelfile);
partTemplate = reshape(rectModel.belta, 5, 5, 31, 9);


max_scale = config.max_scale;
min_scale = config.min_scale;
for i = 1:6
    if ~hasBuff
        img = double(sepSceneSix(i).img .* 255);
        feature = feature_pyramid( img, max_scale, min_scale );
        all_hog{i} = feature;
    else
        feature = all_hog{i};
    end
    for pid = 1:9
        for sid = 1:length(feature.feat)
            match{i, pid, sid} = fconvblas(feature.feat{sid}, {partTemplate(:,:,:,pid)}, 1, 1);
            match{i, pid, sid} = match{i, pid, sid}{1};
        end
    end
end

%% checkrect
% hyps = allhyps(checkrect);
ptIDofSurf = [8 5 1 4; 6 7 3 2; 5 6 2 1; 7 8 4 3; 1 2 3 4; 8 7 6 5];
PROJECTVIEW = [2 1 4 3 6 5];
cuboid_rect_score = cell(length(allhyps),1);
for hid = find(checkrect')
%     fprintf('HID: %d\n', hid);
    hyps = allhyps{hid};
    score = repmat(struct('surfscore',[],'surfID',[]), length(hyps), 1);
    for oid = 1:length(hyps)
%         fprintf('OID: %d\n', oid);
        p1 = hyps(oid).align(1,:);
        p7 = hyps(oid).align(7,:);
        if hyps(oid).type == 1
            tmpstd = std(hyps(oid).align, 1);
            [~,id] = min(tmpstd);
            priorvisibility = false(6,1);
            priorvisibility([id*2-1 id*2]) = true;
        else
            priorvisibility = true(6,1);
        end
        visibility = find(sum(vp .* [p1;p7;p1;p7;p1;p7], 2)<0 & priorvisibility);
        surfscore = -100*ones(length(visibility),1);
        for i = 1:length(visibility)
            scalescore = -100*ones(length(feature.feat),1);
            sid = PROJECTVIEW(visibility(i));
            pid = ptIDofSurf(visibility(i),:);
            [out2D, valid, ~] = projectPoint2SeparateView( hyps(oid).align(pid,:), vp(sid,:), fov, cutSize );
            minscalex1 = 2*8/(out2D(1,1)-1);
            minscalex2 = 3*8/(cutSize - out2D(3,1));
            minscaley1 = 2*8/(out2D(1,2)-1);
            minscaley2 = 3*8/(cutSize - out2D(3,2));
            valid_scale_ID = find(feature.scale>max([minscalex1 minscalex2 minscaley1 minscaley2]));
            if min([minscalex1 minscalex2 minscaley1 minscaley2])<0
                valid_scale_ID = [];
            end

            for j = valid_scale_ID%1:length(feature.feat)
                [xIDmin, yIDmin] = getGridID(out2D(1,1), out2D(1,2), feature.scale(j));
                [xIDmax, yIDmax] = getGridID(out2D(3,1)+1, out2D(3,2)+1, feature.scale(j));
                xID1 = xIDmin;  yID1 = yIDmin;
                xID2 = round((xIDmin+xIDmax)/2);    yID2 = yIDmin;
                xID3 = xIDmax;    yID3 = yIDmin;
                xID4 = xIDmin;  yID4 = round((yIDmin+yIDmax)/2);
                xID5 = round((xIDmin+xIDmax)/2);    yID5 = round((yIDmin+yIDmax)/2);
                xID6 = xIDmax;    yID6 = round((yIDmin+yIDmax)/2);
                xID7 = xIDmin;  yID7 = yIDmax;
                xID8 = round((xIDmin+xIDmax)/2);    yID8 = yIDmax;
                xID9 = xIDmax;    yID9 = yIDmax;
                
                scalescore(j) = match{sid,1,j}(yID1,xID1) + match{sid,2,j}(yID2,xID2) + match{sid,3,j}(yID3,xID3) ...
                              + match{sid,4,j}(yID4,xID4) + match{sid,5,j}(yID5,xID5) + match{sid,6,j}(yID6,xID6) ...
                              + match{sid,7,j}(yID7,xID7) + match{sid,8,j}(yID8,xID8) + match{sid,9,j}(yID9,xID9);
            end
            surfscore(i) = max(scalescore) - rectModel.b;
        end
        score(oid).surfscore = surfscore;
        score(oid).surfID = PROJECTVIEW(visibility);
    end
    cuboid_rect_score{hid} = score;
end


clear rectModel

%% compute multiple layer of segmentation
rotImg  = im2double(rotImg);
thresholds = [200 500 800 1200 2000 3000];
if ~hasBuff
    for i = 1:length(thresholds)
        all_seg(i).labelMap = gbPanoSegment( im2uint8(rotImg), 0.5, thresholds(i), 50);
        all_seg(i).thresh = thresholds(i);
    end
end

%% checkseg
% load('./rectangleDetector/segmentation/uniformvector_lvl6.mat');
[ coor, tri ] = getUniformVector( 6 );
[sH, sW, ~] = size(all_seg(1).labelMap);
coor2D = uv2coords(xyz2uvN(coor), sW, sH);
coor2D_indi = sub2ind([sH sW], coor2D(:,2), coor2D(:,1));


segment_rect_score = cell(length(allhyps),1);
for hid = find(checkseg')
    hyps = allhyps{hid};
    vector_num = size(coor, 1);
    indicator = false(vector_num, length(hyps));
    for i = 2:length(hyps)
        if hyps(i).type==1
            p = hyps(i).out_points_w([1 2 3 4],:);
        else
            p = hyps(i).align;
        end
        p = p./repmat(sqrt(sum(p.^2,2)),1,3);
        vp = sum(p,1); vp = vp./norm(vp); %vp(3) = 0; 
        [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
        if any(~valid)
            vp = sum(p,1); vp(1:2) = 0; vp = vp./norm(vp); 
            [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
            if any(~valid)
                vp = sum(p,1); vp(3) = 0; vp = vp./norm(vp); 
                [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
            end
        end
        K = convhull( out2D(:,1), out2D(:,2));
        [ inside, ~, ~ ] = insideCone( p(K(end:-1:2),:), coor, 0 );
        if ~any(inside)
            fprintf('%d:%d\n', hid, i);
        end
%         assert(any(inside));
        indicator(:,i) = inside;
    end
    
    seg_score = zeros( length(all_seg), length(hyps));
    for sid = 1:length(all_seg)
        label = all_seg(sid).labelMap(coor2D_indi);
        MINLAB = min(label);
        MAXLAB = max(label);
        N = hist(label, MINLAB:MAXLAB);
        for j = 1:length(hyps)
            lab = label(indicator(:,j));
            if isempty(lab)
                continue;
            end
            M = hist(lab, MINLAB:MAXLAB);
            I = M>0.05*sum(M);
            seg_score(sid,j) = sum(M(I))./sum(N(I));
        end
    end
    max_seg_score = max(seg_score, [], 1);
    segment_rect_score{hid} = max_seg_score;
end

if ~hasBuff
    save(bufName, 'all_hog', 'all_seg');
end

end

function [xID, yID] = getGridID(x, y, s)
% xID = round( (x-1) * s / 8 ) + 1;
% yID = round( (y-1) * s / 8 ) + 1;
xID = ceil( (x-1) * s / 8 ) - 2;
yID = ceil( (y-1) * s / 8 ) - 2;
end